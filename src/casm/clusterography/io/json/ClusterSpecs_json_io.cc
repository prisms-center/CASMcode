#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser.hh"
#include "casm/clusterography/io/json/ClusterSpecs_json_io.hh"
#include "casm/clusterography/io/json/IntegralCluster_json_io.hh"
//#include "casm/symmetry/OrbitGeneration_impl.hh"

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"

namespace CASM {

  namespace ClusterSpecs_json_io_impl {

    /// Read vector from 'bspecs' JSON
    ///
    /// \returns std::vector<double> giving <attrname> for clusters in branch 0, 1, etc.
    ///
    template<typename SpecsType>
    std::vector<double> parse_orbit_branch_specs_attr(InputParser<SpecsType> &parser, const std::string &attrname) {

      // "orbit_branch_specs": {
      //   "2": {
      //     "max_length": 6.0
      //   },
      //   "3": {
      //     "max_length": 6.0
      //   }
      // },
      auto specs_it = parser.self.find("orbit_branch_specs");
      if(specs_it == parser.self.end()) {
        parser.error.insert("Error: missing required option 'orbit_branch_specs'.");
        return std::vector<double> {};
      }

      std::vector<double> result;
      auto update_result = [&](int branch, double value) {
        while(branch >= result.size()) {
          result.push_back(0.0);
        }
        result[branch] = value;
      };

      for(auto it = specs_it->begin(); it != specs_it->end(); ++it) {
        auto attr_it = it->find(attrname);
        if(attr_it != it->end()) {
          update_result(std::stoi(it.name()), attr_it->template get<double>());
        }
      }
      return result;
    }

    void write_group_indices(
      const SymGroup &grp,
      jsonParser &json) {

      json.put_array();
      for(const auto &op : grp) {
        json.push_back(op.index());
      }
    }

    std::vector<SymOp>::const_iterator find_by_master_group_index(
      const SymOp &op,
      const SymGroup &super_group) {
      auto it = super_group.begin();
      auto end = super_group.end();
      for(; it != end; ++it) {
        if(it->index() == op.index()) {
          return it;
        }
      }
      return end;
    }

    std::unique_ptr<SymGroup> group_from_indices(
      const SymGroup &super_group,
      const std::vector<Index> &indices) {

      auto result = notstd::make_unique<SymGroup>();
      result->set_lattice(super_group.lattice());
      const MasterSymGroup &master_group = super_group.master_group();
      for(const auto &i : indices) {
        const auto &op = master_group[i];

        auto it = find_by_master_group_index(op, super_group);
        if(it == super_group.end()) {
          continue;
        }
        result->push_back(*it);
      }
      result->sort();
      if((*result)[0].index() != 0) {
        throw std::runtime_error("Error in group_from_indices: First element is not identity.");
      }
      if(!result->is_group(super_group.lattice().tol())) {
        throw std::runtime_error("Error in local_group_from_indices: Not a group.");
      }
      return result;
    }

    std::unique_ptr<SymGroup> local_group_from_indices(
      const SymGroup &super_group,
      const std::vector<Index> &indices,
      const IntegralCluster &phenomenal,
      const PrimPeriodicSymCompare<IntegralCluster> &sym_compare) { // TODO: check this?

      IntegralCluster e {sym_compare.prepare(phenomenal)};

      const MasterSymGroup &master_group = super_group.master_group();

      auto result = notstd::make_unique<SymGroup>();
      result->set_lattice(master_group.lattice());
      for(const auto &i : indices) {

        const auto &op = master_group[i];

        auto it = find_by_master_group_index(op, super_group);
        if(it == super_group.end()) {
          continue;
        }

        // get translation & check that phenomenal cluster sites remain invariant
        // - note: this is a minimum requirement check that provides the translation, it does not
        //   guarantee the user input op is correct for a particular context (i.e. DiffTrans)
        if(sym_compare.equal(e, sym_compare.prepare(sym_compare.copy_apply(op, e)))) {
          result->push_back(sym_compare.spatial_transform()*op);
        }
        else {
          throw std::runtime_error("Error in local_group_from_indices: Phenomenal cluster sites are not invariant.");
        }
      }
      result->sort();
      if((*result)[0].index() != 0) {
        throw std::runtime_error("Error in local_group_from_indices: First element is not identity.");
      }
      if(!result->is_group(sym_compare.tol())) {
        throw std::runtime_error("Error in local_group_from_indices: Not a group.");
      }
      return result;
    }

    std::unique_ptr<SymGroup> parse_generating_group(
      InputParser<PeriodicMaxLengthClusterSpecs> &parser,
      const std::shared_ptr<const Structure> &shared_prim,
      const SymGroup &super_group) {

      auto generating_group_indices = parser.optional<std::vector<Index>>("generating_group");
      if(generating_group_indices) {
        try {
          return group_from_indices(super_group, *generating_group_indices);
        }
        catch(std::exception &e) {
          parser.error.insert(std::string("Error parsing generating_group:") + e.what());
          return std::unique_ptr<SymGroup> {};
        }
      }
      else {
        return notstd::clone(super_group);
      }
    }

    std::unique_ptr<SymGroup> parse_local_generating_group(
      InputParser<LocalMaxLengthClusterSpecs> &parser,
      const std::shared_ptr<const Structure> &shared_prim,
      const IntegralCluster &phenomenal,
      const SymGroup &super_group) {

      std::vector<Index> generating_group_indices;
      parser.require(generating_group_indices, "generating_group");
      if(!parser.valid()) {
        return std::unique_ptr<SymGroup>();
      }
      try {
        typedef PrimPeriodicSymCompare<IntegralCluster> SymCompareType;
        SymCompareType sym_compare {shared_prim, shared_prim->lattice().tol()}; // TODO: check this?
        return local_group_from_indices(
                 super_group,
                 generating_group_indices,
                 phenomenal,
                 sym_compare);
      }
      catch(std::exception &e) {
        parser.error.insert(std::string("Error parsing generating_group:") + e.what());
      }
      return std::unique_ptr<SymGroup>();
    }
  }

  /// Parse PeriodicMaxLengthClusterSpecs from JSON
  ///
  /// Expects:
  /// \code
  /// {
  ///   "orbit_branch_specs": { // required, start at "2"
  ///     "2": {"max_length":<number>, "cutoff_radius":<number>},
  ///     "3": {"max_length":<number>, "cutoff_radius":<number>},
  ///     ...},
  ///   "orbit_specs": [ // optional
  ///      ... prototype periodic clusters...
  ///   ]
  /// }
  /// \endcode
  ///
  void parse(
    InputParser<PeriodicMaxLengthClusterSpecs> &parser,
    const std::shared_ptr<const Structure> &shared_prim,
    const SymGroup &super_group) {
    parser.value = notstd::make_unique<PeriodicMaxLengthClusterSpecs>(shared_prim);
    auto &cspecs = *parser.value;

    using namespace ClusterSpecs_json_io_impl;
    cspecs.max_length = parse_orbit_branch_specs_attr(parser, "max_length");
    auto generating_group = parse_generating_group(parser, shared_prim, super_group);
    if(generating_group) {
      cspecs.generating_group = std::move(generating_group);
    }
    parser.subparse_if(cspecs.custom_generators, "orbit_specs", shared_prim);
    cspecs.site_filter = alloy_sites_filter; // TODO: update for all dof
  }

  /// Parse LocalMaxLengthClusterSpecs from JSON
  ///
  /// \param parser InputParser stores resulting value, errors, and warnings
  /// \param shared_prim Prim structure
  /// \param super_group Super group for generating sub groups. When reading the generating_group
  ///        from JSON, only members of this group will be retained.
  ///
  /// Expects:
  /// \code
  /// {
  ///   "generating_group": [0, ...], // required, array of symop indices for generating group
  ///   "phenomenal": { // IntegralCluster as JSON
  ///     "coordinate_mode" : <COORD_TYPE>, // ("FRAC", "CART", "INT" (default)) (optional)
  ///     "sites" : [
  ///       [b, i, j, k], // of type matching coordinate_mode
  ///       ...
  ///     ]
  ///   },
  ///   "orbit_branch_specs": { // required, start at "1"
  ///     "1": {"max_length":<number>, "cutoff_radius":<number>},
  ///     "2": {"max_length":<number>, "cutoff_radius":<number>},
  ///     ...},
  ///   "orbit_specs": [ // optional
  ///      ... prototype local clusters...
  ///   ]
  /// }
  /// \endcode
  ///
  ///
  void parse(
    InputParser<LocalMaxLengthClusterSpecs> &parser,
    const std::shared_ptr<const Structure> &shared_prim,
    const SymGroup &super_group) {

    using namespace ClusterSpecs_json_io_impl;

    IntegralCluster phenomenal {*shared_prim};
    parser.require(phenomenal, "phenomenal", shared_prim->lattice().tol());

    auto local_group = parse_local_generating_group(
                         parser,
                         shared_prim,
                         phenomenal,
                         super_group);

    parser.value = notstd::make_unique<LocalMaxLengthClusterSpecs>(
                     shared_prim,
                     std::move(local_group),
                     phenomenal);

    auto &cspecs = *parser.value;

    cspecs.max_length = parse_orbit_branch_specs_attr(parser, "max_length");
    cspecs.cutoff_radius = parse_orbit_branch_specs_attr(parser, "cutoff_radius");
    parser.subparse_if(cspecs.custom_generators, "orbit_specs", shared_prim);
    cspecs.site_filter = alloy_sites_filter; // TODO: update for all dof
  }


  /// \brief Parse PeriodicMaxLengthClusterSpecs or LocalMaxLengthClusterSpecs from JSON & validate
  ///
  /// This is used by PrimClex to read any user's bspecs.json file. Update this if you want to add
  /// more recognized cluster specs.
  ///
  /// This overload uses `shared_prim->factor_group()` for the generating group.
  ///
  /// Expects:
  /// \code
  /// {
  ///   "method": <method name>, // one of "periodic_max_length", "local_max_length"
  ///   "params": { ... depends on method ... }
  /// }
  /// \endcode
  ///
  /// For method == "periodic_max_length", "params" is:
  /// \code
  /// {
  ///   "orbit_branch_specs": { // required, start at "2"
  ///     "2": {"max_length":<number>, "cutoff_radius":<number>},
  ///     "3": {"max_length":<number>, "cutoff_radius":<number>},
  ///     ...},
  ///   "orbit_specs": [ // optional
  ///      ... prototype periodic clusters...
  ///   ]
  /// }
  /// \endcode
  ///
  ///
  /// For method == "local_max_length", "params" is:
  /// \code
  /// {
  ///   "diff_trans": <PrimPeriodicDiffTransOrbit name>, // required, only "diff_trans" for now
  ///   "orbit_branch_specs": { // required, start at "1"
  ///     "1": {"max_length":<number>, "cutoff_radius":<number>},
  ///     "2": {"max_length":<number>, "cutoff_radius":<number>},
  ///     ...},
  ///   "orbit_specs": [ // optional
  ///      ... prototype local clusters...
  ///   ]
  /// }
  /// \endcode
  ///
  void parse(
    InputParser<ClusterSpecs> &parser,
    const std::shared_ptr<const Structure> &shared_prim) {
    parse(parser, shared_prim, shared_prim->factor_group());
  }

  /// \brief Parse PeriodicMaxLengthClusterSpecs or LocalMaxLengthClusterSpecs from JSON & validate with specified generating_group
  ///
  /// This is used by PrimClex to read any user's bspecs.json file. Update this if you want to add
  /// more recognized cluster specs.
  ///
  /// Expects:
  /// \code
  /// {
  ///   "method": <method name>, // one of "periodic_max_length", "local_max_length"
  ///   "params": { ... depends on method ... }
  /// }
  /// \endcode
  ///
  /// For method == "periodic_max_length", "params" is:
  /// \code
  /// {
  ///   "orbit_branch_specs": { // required, start at "2"
  ///     "2": {"max_length":<number>, "cutoff_radius":<number>},
  ///     "3": {"max_length":<number>, "cutoff_radius":<number>},
  ///     ...},
  ///   "orbit_specs": [ // optional
  ///      ... prototype periodic clusters...
  ///   ]
  /// }
  /// \endcode
  ///
  ///
  /// For method == "local_max_length", "params" is:
  /// \code
  /// {
  ///   "diff_trans": <PrimPeriodicDiffTransOrbit name>, // required, only "diff_trans" for now
  ///   "orbit_branch_specs": { // required, start at "1"
  ///     "1": {"max_length":<number>, "cutoff_radius":<number>},
  ///     "2": {"max_length":<number>, "cutoff_radius":<number>},
  ///     ...},
  ///   "orbit_specs": [ // optional
  ///      ... prototype local clusters...
  ///   ]
  /// }
  /// \endcode
  ///
  void parse(
    InputParser<ClusterSpecs> &parser,
    const std::shared_ptr<const Structure> &shared_prim,
    const SymGroup &super_group) {

    std::string method;
    parser.require(method, "method");
    if(method.empty()) {
      return;
    }

    // ** this could be a dictionary lookup **
    if(method == "periodic_max_length") {
      parser.subparse_as<PeriodicMaxLengthClusterSpecs>("params", shared_prim, super_group);
    }
    else if(method == "local_max_length") {
      parser.subparse_as<LocalMaxLengthClusterSpecs>("params", shared_prim, super_group);
    }
    else {
      parser.error.insert("Error: unknown cluster_specs method '" + method + "'.");
    }
  }


  namespace ClusterSpecs_json_io_impl {

    jsonParser &orbit_branch_specs_attr_to_json(std::vector<double> attr, const std::string &attrname, jsonParser &json) {
      jsonParser &j = json["orbit_branch_specs"];
      for(int b = 0; b < attr.size(); ++b) {
        j[std::to_string(b)][attrname] = attr[b];
      }
      return json;
    }
  }

  jsonParser &to_json(
    const PeriodicMaxLengthClusterSpecs &cspecs,
    jsonParser &json) {

    json["method"] = cspecs.name();
    json["params"] = jsonParser::object();

    using namespace ClusterSpecs_json_io_impl;
    orbit_branch_specs_attr_to_json(cspecs.max_length, "max_length", json["params"]);
    write_group_indices(*cspecs.generating_group, json["params"]["generating_group"]);
    return json;
  }


  /// \brief Write LocalMaxLengthClusterSpecs to JSON
  jsonParser &to_json(
    const LocalMaxLengthClusterSpecs &cspecs,
    jsonParser &json) {

    json["method"] = cspecs.name();
    json["params"] = jsonParser::object();

    // select based on method
    using namespace ClusterSpecs_json_io_impl;
    orbit_branch_specs_attr_to_json(cspecs.max_length, "max_length", json["params"]);
    orbit_branch_specs_attr_to_json(cspecs.cutoff_radius, "cutoff_radius", json["params"]);
    write_group_indices(*cspecs.generating_group, json["params"]["generating_group"]);
    to_json(cspecs.phenomenal, json["params"]["phenomenal"]);
    return json;
  }

  /// \brief Write PeriodicMaxLengthClusterSpecs or LocalMaxLengthClusterSpecs to JSON
  jsonParser &to_json(
    const ClusterSpecs &cspecs,
    jsonParser &json) {
    if(cspecs.name() == PeriodicMaxLengthClusterSpecs::method_name) {
      return to_json(static_cast<const PeriodicMaxLengthClusterSpecs &>(cspecs), json);
    }
    else if(cspecs.name() == LocalMaxLengthClusterSpecs::method_name) {
      return to_json(static_cast<const LocalMaxLengthClusterSpecs &>(cspecs), json);
    }
    else {
      throw std::runtime_error(std::string("Error converting ClusterSpecs to JSON:") + " cannot convert method: '" + cspecs.name() + "'");
    }
  }

}
