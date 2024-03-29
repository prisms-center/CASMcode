#include "casm/clusterography/io/json/ClusterSpecs_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/clusterography/SupercellClusterOrbits.hh"
#include "casm/clusterography/io/json/ClusterOrbits_json_io.hh"
#include "casm/clusterography/io/json/IntegralCluster_json_io.hh"

namespace CASM {

namespace ClusterSpecs_json_io_impl {

/// Read vector from 'bspecs' JSON
///
/// \returns std::vector<double> giving <attrname> value branch 0, 1, etc.
///
/// Expects `parser.self` in format:
/// \code
/// "orbit_branch_specs": {
///   // includes "0": {<attrname>: 0.0} always, by default
///   // includes "b": {<attrname>: 0.0} if any greater branch is present
///   "2": { <attrname>: 8.0 },
///   "3": { <attrname>: 6.0 }
/// }
/// \endcode
///
template <typename SpecsType>
std::vector<double> parse_orbit_branch_specs_attr(
    InputParser<SpecsType> &parser, const std::string &attrname) {
  auto specs_it = parser.self.find("orbit_branch_specs");
  if (specs_it == parser.self.end() || specs_it->is_null()) {
    return std::vector<double>({0.});
  }
  std::vector<double> result({0.});

  if (!specs_it->is_obj()) {
    fs::path path{"orbit_branch_specs"};
    std::stringstream msg;
    msg << "Expected a JSON object";
    parser.insert_error("orbit_branch_specs", msg.str());
    return result;
  }

  auto update_result = [&](int branch, double value) {
    while (branch >= result.size()) {
      result.push_back(0.0);
    }
    result[branch] = value;
  };

  for (auto it = specs_it->begin(); it != specs_it->end(); ++it) {
    double value = 0.;  // default is 0. if no attrname found
    if (!it->is_obj()) {
      fs::path path{"orbit_branch_specs"};
      std::stringstream msg;
      msg << "Expected JSON object containing '" << attrname
          << "' (double, default=0.).";
      parser.insert_error(path / it.name(), msg.str());
      continue;
    }
    auto attr_it = it->find(attrname);
    if (attr_it != it->end()) {
      value = attr_it->template get<double>();
    }
    update_result(std::stoi(it.name()), value);
  }
  return result;
}

void write_group_indices(const SymGroup &grp, jsonParser &json) {
  json.put_array();
  for (const auto &op : grp) {
    json.push_back(op.index());
  }
}

void write_group_indices(std::vector<PermuteIterator> const &grp,
                         jsonParser &json) {
  std::set<Index> prim_factor_group_indices;
  for (auto const &permute_it : grp) {
    prim_factor_group_indices.insert(permute_it.prim_factor_group_index());
  }

  json = prim_factor_group_indices;
}

std::vector<SymOp>::const_iterator find_by_master_group_index(
    const SymOp &op, const SymGroup &super_group) {
  auto it = super_group.begin();
  auto end = super_group.end();
  for (; it != end; ++it) {
    if (it->index() == op.index()) {
      return it;
    }
  }
  return end;
}

std::unique_ptr<SymGroup> group_from_indices(const SymGroup &super_group,
                                             const std::set<Index> &indices) {
  auto result = notstd::make_unique<SymGroup>();
  result->set_lattice(super_group.lattice());
  for (SymOp const &op : super_group) {
    if (indices.count(op.index())) {
      result->push_back(op);
    }
  }
  result->sort();
  if ((*result)[0].index() != 0) {
    throw std::runtime_error(
        "Error in group_from_indices: First element is not identity.");
  }
  if (!result->is_group(super_group.lattice().tol())) {
    throw std::runtime_error("Error in local_group_from_indices: Not a group.");
  }
  return result;
}

std::unique_ptr<SymGroup> local_group_from_indices(
    const SymGroup &super_group,
    const std::set<Index> &phenomenal_group_indices,
    const IntegralCluster &phenomenal,
    std::shared_ptr<const Structure> const &shared_prim) {
  PrimPeriodicSymCompare<IntegralCluster> sym_compare{
      shared_prim, shared_prim->lattice().tol()};
  IntegralCluster e{sym_compare.prepare(phenomenal)};
  auto result = notstd::make_unique<SymGroup>();
  result->set_lattice(super_group.lattice());
  for (SymOp const &op : super_group) {
    if (phenomenal_group_indices.count(op.index())) {
      // get translation & check that phenomenal cluster sites remain invariant
      // - note: this is a minimum requirement check that provides the
      // translation, it does not
      //   guarantee the user input indices are correct for a particular context
      //   (i.e. DiffTrans)
      if (sym_compare.equal(
              e, sym_compare.prepare(sym_compare.copy_apply(op, e)))) {
        result->push_back(sym_compare.spatial_transform() * op);
      } else {
        throw std::runtime_error(
            "Error in local_group_from_indices: Phenomenal cluster sites are "
            "not invariant.");
      }
    }
  }
  result->sort();
  if ((*result)[0].index() != 0) {
    throw std::runtime_error(
        "Error in local_group_from_indices: First element is not identity.");
  }
  if (!result->is_group(sym_compare.tol())) {
    throw std::runtime_error("Error in local_group_from_indices: Not a group.");
  }
  return result;
}

std::unique_ptr<SymGroup> parse_generating_group(
    InputParser<PeriodicMaxLengthClusterSpecs> &parser,
    const std::shared_ptr<const Structure> &shared_prim,
    const SymGroup &super_group) {
  auto generating_group_indices =
      parser.optional<std::set<Index>>("generating_group");
  if (generating_group_indices) {
    try {
      return group_from_indices(super_group, *generating_group_indices);
    } catch (std::exception &e) {
      parser.error.insert(std::string("Error parsing generating_group:") +
                          e.what());
      return std::unique_ptr<SymGroup>{};
    }
  } else {
    return notstd::clone(super_group);
  }
}

std::unique_ptr<SymGroup> parse_local_generating_group(
    InputParser<LocalMaxLengthClusterSpecs> &parser,
    const std::shared_ptr<const Structure> &shared_prim,
    const IntegralCluster &phenomenal, const SymGroup &super_group) {
  std::set<Index> generating_group_indices;
  parser.require(generating_group_indices, "generating_group");
  if (!parser.valid()) {
    return std::unique_ptr<SymGroup>();
  }
  try {
    return local_group_from_indices(super_group, generating_group_indices,
                                    phenomenal, shared_prim);
  } catch (std::exception &e) {
    parser.error.insert(std::string("Error parsing generating_group:") +
                        e.what());
  }
  return std::unique_ptr<SymGroup>();
}
}  // namespace ClusterSpecs_json_io_impl

/// Parse PeriodicMaxLengthClusterSpecs from JSON
///
/// Expects:
/// \code
/// {
///   "orbit_branch_specs": { // optional
///     "2": {"max_length":<number>, "cutoff_radius":<number>},
///     "3": {"max_length":<number>, "cutoff_radius":<number>},
///     ...},
///   "orbit_specs": [ // optional
///      ... prototype periodic clusters...
///   ]
/// }
/// \endcode
///
/// Notes:
/// - The null cluster orbit ("0") is always included with max_length value 0.0
///   (the value itself does not matter).
/// - The point cluster orbit ("1") is included with max_length value 0.0 (the
///   value itself does not matter), if "1" is present with any value, or if
///   "2" or greater is included.
///
void parse(InputParser<PeriodicMaxLengthClusterSpecs> &parser,
           const std::shared_ptr<const Structure> &shared_prim,
           const SymGroup &super_group) {
  using namespace ClusterSpecs_json_io_impl;

  // parse generating group
  auto generating_group_ptr =
      parse_generating_group(parser, shared_prim, super_group);
  if (!generating_group_ptr) {
    generating_group_ptr = notstd::clone(super_group);
  }

  // parse max length
  auto max_length = parse_orbit_branch_specs_attr(parser, "max_length");

  // parse custom generators ("orbit_specs")
  std::vector<IntegralClusterOrbitGenerator> default_custom_generators{};
  auto custom_generators_parser =
      parser.subparse_else<std::vector<IntegralClusterOrbitGenerator>>(
          "orbit_specs", default_custom_generators, *shared_prim);

  if (!parser.valid()) {
    return;
  }
  parser.value = notstd::make_unique<PeriodicMaxLengthClusterSpecs>(
      shared_prim, *generating_group_ptr, dof_sites_filter(), max_length,
      *custom_generators_parser->value);
}

/// Parse LocalMaxLengthClusterSpecs from JSON
///
/// \param parser InputParser stores resulting value, errors, and warnings
/// \param shared_prim Prim structure
/// \param super_group Super group for generating sub groups. When reading the
/// generating_group
///        from JSON, only members of this group will be retained.
///
/// Expects:
/// \code
/// {
///   "generating_group": [0, ...], // required, array of symop indices for
///   generating group "phenomenal": { // IntegralCluster as JSON
///     "coordinate_mode" : <COORD_TYPE>, // ("FRAC", "CART", "INT" (default))
///     (optional) "sites" : [
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
void parse(InputParser<LocalMaxLengthClusterSpecs> &parser,
           const std::shared_ptr<const Structure> &shared_prim,
           const SymGroup &super_group) {
  using namespace ClusterSpecs_json_io_impl;

  // parse phenomenal
  auto phenomenal_subparser_ptr =
      parser.subparse<IntegralCluster>("phenomenal", *shared_prim);
  if (!phenomenal_subparser_ptr->valid()) {
    return;
  }
  IntegralCluster phenomenal = *phenomenal_subparser_ptr->value;

  // parse generating group
  auto generating_group_ptr = parse_local_generating_group(
      parser, shared_prim, phenomenal, super_group);

  // parse max_length and cutoff_radius
  auto max_length = parse_orbit_branch_specs_attr(parser, "max_length");
  auto cutoff_radius = parse_orbit_branch_specs_attr(parser, "cutoff_radius");

  // parse custom generators ("orbit_specs")
  std::vector<IntegralClusterOrbitGenerator> default_custom_generators{};
  auto custom_generators_parser =
      parser.subparse_else<std::vector<IntegralClusterOrbitGenerator>>(
          "orbit_specs", default_custom_generators, *shared_prim);

  // TODO: include option in JSON?
  bool include_phenomenal_sites = false;

  if (!parser.valid()) {
    return;
  }
  parser.value = notstd::make_unique<LocalMaxLengthClusterSpecs>(
      shared_prim, *generating_group_ptr, phenomenal, dof_sites_filter(),
      max_length, cutoff_radius, include_phenomenal_sites,
      *custom_generators_parser->value);
}

/// \brief Parse PeriodicMaxLengthClusterSpecs or LocalMaxLengthClusterSpecs
/// from JSON & validate
///
/// This overload is equivalent to: `parse(parser, shared_prim,
/// shared_prim->factor_group())`.
///
void parse(InputParser<ClusterSpecs> &parser,
           const std::shared_ptr<const Structure> &shared_prim) {
  parse(parser, shared_prim, shared_prim->factor_group());
}

/// \brief Parse any supported type of ClusterSpecs from JSON & validate with
/// specified generating_group
///
/// \param parser (InputParser<ClusterSpecs>) Holds resulting value and error or
/// warning messages \param shared_prim (const std::shared_ptr<const Structure>
/// &) Prim structure \param super_group (const SymGroup &) Used to generate
/// orbits.
/// \param transformation_matrix_to_super (Eigen::Matrix3l const *) Optional
/// super lattice
///        transformation matrix for contexts when "within_scel_max_length" is
///        supported.
///
/// This is used by PrimClex to read any user's bspecs.json file. Update this if
/// you want to add more recognized cluster specs.
///
/// Notes on the naming conventions:
/// - An "orbit" is the set of all equivalent objects under some generating
///   symmetry group. It is a general term, and in different contexts within
///   CASM could be used for orbits of clusters, orbits of functions, or orbits
///   of other objects.
/// - There are infinitely many orbits of 2-point clusters, corresponding to 1st
///   nearest neighbor, 2nd nearest neighbor, 3rd nearest neighbor, etc.
///   clusters. The same applies to 3-point, 4-point, etc. clusters.
/// - The "orbit tree" of all clusters is made up of "orbit branches" where each
///   "branch" of the "orbit tree" is made up of all clusters with the same
///   number of sites. For example, the n=2 "orbit branch" consists of all
///   orbits of 2-point clusters, the n=3 "orbit branch" consists of all orbits
///   of 3-point clusters, etc.
///
/// Notes on cluster generation methods:
/// - Clusters are generated recursively, meaning 2-point clusters are generated
///   by attempting to add sites to 1-point clusters, 3-point clusters are
///   generated by attempting to add sites to 2-point clusters, etc. This means
///   that some (n+1)-point clusters that would otherwise be included may not
///   be if the (n+1)-point cluster `max_length` is greater than the n-point
///   cluster `max_length`. (The max_length value does not matter for null and
///   1-point cluster orbits, so this applies to n >= 2).
///
/// Expected JSON format:
///
///     method: string (required)
///         Specify which cluster orbit generating method will be used. One of:
///         - "periodic_max_length": Clusters differing by a lattice translation
///         are considered
///           equivalent. Cluster generation is truncated by specifying the
///           maximum distance between sites in a cluster for 2-point, 3-point,
///           etc. clusters. The point clusters comprising the asymmetric unit
///           of the prim structure are always included. This is the appropriate
///           method for cluster expansions of crystal properties.
///         - "local_max_length": Clusters differing by a lattice translation
///         are considered
///           distinct. Cluster generation is truncated by specifying for
///           2-point, 3-point, etc. clusters the radius of sites from
///           "phenomenal" cluster include, and the maximum distance between
///           sites in a cluster. This is the appropriate method for cluster
///           expansions of properties associated with the "phenomenal" cluster.
///     params: object (required)
///         Specifies parameters for the method selected by `method`. Options
///         depend on the `method` chosen:
///
///         For method=="periodic_max_length":
///             orbit_branch_specs: object (required)
///                 Cluster generation is truncated by specifying the maximum
///                 distance between sites in a cluster for each orbit branch
///                 (i.e. 2-point, 3-point, etc. clusters). The null (0-point)
///                 clusters are always included. The 1-point clusters
///                 comprising the asymmetric unit of the prim structure are
///                 always included if 2-point clusters or greater are included
///                 or if "1" is present (with any value).
///
///
///                 Example 1: Generate up to 3-point clusters.
///                 - null cluster orbit is always included
///                 - 1-point cluster orbits are included because 2-point
///                   cluster orbits are included
///                 - max_length value for point branch 3 is >= max_length
///                   value for branch 2
///
///                     "orbit_branch_specs": {
///                         "2": { "max_length": 10.0 },
///                         "3": { "max_length": 8.0 }
///                     }
///
///                 Example 2: Generate null cluster orbit only
///                 - null cluster orbit is always included
///                   ("orbit_branch_specs" itself need not be included)
///
///                     "orbit_branch_specs": {}
///
///                 Example 3: Generate null and 1-point cluster orbits only
///                 - null cluster orbit is always included
///                 - only the presence of "1" needed to included 1-point
///                   cluster orbits, the max_length value is not necessary for
///                   1-point clusters and ignored if present
///
///                     "orbit_branch_specs": {
///                         "1": {}
///                     }
///
///             orbit_specs: array (optional)
///                 An array of clusters which are used to generate and include
///                 orbits of clusters whether or not they meet the `max_length`
///                 truncation criteria. See the cluster input format below.
///
///         For method=="local_max_length":
///             phenomenal: object (required)
///                 The "phenomenal" cluster about which local clusters are
///                 generated. See the cluster input format below.
///             generating_group: array of int (required)
///                 An array of symop indices into the prim structure factor
///                 group specifying the invariant group of the "phenomenal"
///                 cluster which should be used to for generating local cluster
///                 orbits. In some contexts, the relevant symmetry is lower
///                 than that determined from the phenomenal cluster sites
///                 alone.
///             orbit_branch_specs: object (required)
///                 All sites within `cutoff_radius` distance of any site in the
///                 phenomenal cluster are considered candidates for inclusion
///                 in clusters of a particular orbit branch. Cluster generation
///                 is truncated by specifying the maximum distance between
///                 sites in a cluster for each orbit branch. The `max_length`
///                 parameter is not necessary for 1-point clusters and ignored
///                 if present.
///
///                 Example:
///                     "orbit_branch_specs": {
///                         "1": { "cutoff_radius": 6.0 },
///                         "2": { "max_length": 9.0, "cutoff_radius": 6.0 },
///                         "3": { "max_length": 8.0, "cutoff_radius": 6.0 },
///                         ...
///                     }
///
///             orbit_specs: array (optional)
///                 An array of clusters which are used to generate and include
///                 orbits of clusters whether or not they meet the
///                 `cutoff_radius` or `max_length` truncation criteria. See the
///                 cluster input format below.
///
/// Cluster input format for "orbit_specs" and "phenomenal": object
///     coordinate_mode: string (optional, default="Integral")
///         Specifies the coordinate mode used to specify cluster sites. May be:
///         - "Integral": 4-index coordinate [b, i, j, k], where b=sublattice
///         index, and i,j,k
///           are lattice vector indices. Also accepts "INT", "INTEGRAL",
///           "integral".
///         - "Cartesian": 3-index coordinate [x, y, z] giving the site in
///         Cartesian coordinates.
///           Also accepts "CART", "cartesian".
///         - "Direct" or "Fractional": 3-index coordinate [a, b, c], where
///         a,b,c are multiplied
///           by the lattice vectors to give the site coordinate in Cartesian
///           coordinates. Also accepts "FRAC", "fractional", "direct".
///     sites: array of arrays (required)
///         An array of coordinates of sites in the cluster.
///     prototype: this is an allowed alias for "sites"
///     include_subclusters: boolean (optional, default=true)
///         Whether all subclusters of the specified clusters should also be
///         included. This is not relevant for "phenomenal" and ignored if
///         present.
///
///     Example cluster, with "Direct" coordinates:
///         {
///             "coordinate_mode" : "Direct",
///             "sites" : [
///                 [ 0.000000000000, 0.000000000000, 0.000000000000 ],
///                 [ 1.000000000000, 0.000000000000, 0.000000000000 ],
///                 [ 2.000000000000, 0.000000000000, 0.000000000000 ],
///                 [ 3.000000000000, 0.000000000000, 0.000000000000 ]],
///             "include_subclusters" : true
///         }
///
///     Example cluster, with "Integral" coordinates:
///         {
///             "coordinate_mode" : "Integral",
///             "sites" : [
///                 [ 0, 0, 0, 0 ],
///                 [ 0, 1, 0, 0 ],
///                 [ 1, 0, 0, 0 ]],
///             "include_subclusters" : true
///         }
///
void parse(InputParser<ClusterSpecs> &parser,
           const std::shared_ptr<const Structure> &shared_prim,
           const SymGroup &super_group) {
  std::string method;
  parser.require(method, "method");
  if (method.empty()) {
    return;
  }

  // ** this could be a dictionary lookup **
  if (method == "periodic_max_length") {
    auto subparser = parser.subparse<PeriodicMaxLengthClusterSpecs>(
        "params", shared_prim, super_group);
    if (subparser->value) {
      parser.value = std::move(subparser->value);
    }
  } else if (method == "local_max_length") {
    auto subparser = parser.subparse<LocalMaxLengthClusterSpecs>(
        "params", shared_prim, super_group);
    if (subparser->value) {
      parser.value = std::move(subparser->value);
    }
  } else {
    parser.error.insert("Error: unknown cluster_specs method '" + method +
                        "'.");
  }
}

namespace ClusterSpecs_json_io_impl {

jsonParser &orbit_branch_specs_attr_to_json(std::vector<double> attr,
                                            const std::string &attrname,
                                            jsonParser &json) {
  jsonParser &j = json["orbit_branch_specs"];
  for (int b = 0; b < attr.size(); ++b) {
    j[std::to_string(b)][attrname] = attr[b];
  }
  return json;
}
}  // namespace ClusterSpecs_json_io_impl

jsonParser &to_json(const PeriodicMaxLengthClusterSpecs &cspecs,
                    jsonParser &json) {
  json["method"] = cspecs.name();
  jsonParser &json_params = json["params"];

  using namespace ClusterSpecs_json_io_impl;
  orbit_branch_specs_attr_to_json(cspecs.max_length, "max_length", json_params);
  write_group_indices(cspecs.generating_group, json_params["generating_group"]);
  if (cspecs.custom_generators.size()) {
    json_params["orbit_specs"] = cspecs.custom_generators;
  }
  // currently fixed: site_filter=dof_sites_filter()
  return json;
}

/// \brief Write LocalMaxLengthClusterSpecs to JSON
jsonParser &to_json(const LocalMaxLengthClusterSpecs &cspecs,
                    jsonParser &json) {
  json["method"] = cspecs.name();
  jsonParser &json_params = json["params"];

  // select based on method
  using namespace ClusterSpecs_json_io_impl;
  orbit_branch_specs_attr_to_json(cspecs.max_length, "max_length", json_params);
  orbit_branch_specs_attr_to_json(cspecs.cutoff_radius, "cutoff_radius",
                                  json_params);
  write_group_indices(cspecs.generating_group, json_params["generating_group"]);
  to_json(cspecs.phenomenal, json_params["phenomenal"]);
  if (cspecs.custom_generators.size()) {
    json_params["orbit_specs"] = cspecs.custom_generators;
  }
  // currently fixed: site_filter=dof_sites_filter()
  return json;
}

/// \brief Write PeriodicMaxLengthClusterSpecs or LocalMaxLengthClusterSpecs to
/// JSON
jsonParser &to_json(const ClusterSpecs &cspecs, jsonParser &json) {
  if (cspecs.name() == PeriodicMaxLengthClusterSpecs::method_name) {
    return to_json(static_cast<const PeriodicMaxLengthClusterSpecs &>(cspecs),
                   json);
  } else if (cspecs.name() == LocalMaxLengthClusterSpecs::method_name) {
    return to_json(static_cast<const LocalMaxLengthClusterSpecs &>(cspecs),
                   json);
  } else {
    throw std::runtime_error(
        std::string("Error converting ClusterSpecs to JSON:") +
        " cannot convert method: '" + cspecs.name() + "'");
  }
}

}  // namespace CASM
