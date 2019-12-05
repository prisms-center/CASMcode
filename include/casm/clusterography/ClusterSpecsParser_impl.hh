#ifndef CASM_ClusterSpecsParser_impl
#define CASM_ClusterSpecsParser_impl

#include "casm/clusterography/ClusterSpecsParser.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/global/definitions.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/database/Database.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/database/ConfigDatabase.hh"

namespace CASM {

  template<typename RequiredType>
  bool OrbitBranchSpecsParser::require_previous(
    jsonParser::const_iterator it,
    std::string option) {

    auto prev_it = previous(it);
    if(prev_it == self.end()) {
      std::stringstream msg;
      msg << "Error: "
          << "branch '" << branch_to_string(branch_to_int(it) - 1) << "' is required because branch '"
          << it.name() << "' is included.";
      error.insert(msg.str());
      return false;
    }
    return bool(require<RequiredType>(it, option));
  }

  template<typename RequiredType>
  bool OrbitBranchSpecsParser::require_nonincreasing(
    jsonParser::const_iterator it,
    std::string option) {

    if(!require<RequiredType>(it, option) || !require_previous<RequiredType>(it, option)) {
      return false;
    }

    if(it->find(option)->get<double>() > previous(it)->find(option)->get<double>()) {
      std::stringstream msg;
      msg << "Error: "
          << "'" << option << "' increases from branch '" << previous(it).name()
          << "' to branch '" << it.name() << "'";
      error.insert(msg.str());
      return false;
    }
    return true;
  }

  // --- struct LocalOrbitSpecsParser ---

  /// For all custom clusters, insert 'op*prototype' into custom_generators
  ///
  /// \note Use ClusterEquivalenceParser to determine if custom clusters apply
  /// to a given cluster, and determine 'op'.
  ///
  template<typename OrbitType>
  OrbitGenerators<OrbitType> &LocalOrbitSpecsParser::insert_custom_generators(
    const SymOp &op,
    OrbitGenerators<OrbitType> &custom_generators) const {
    for(const auto &val : prototypes) {
      if(val.include_subclusters) {
        insert_subcluster_generators(copy_apply(op, val.cluster), custom_generators, primclex.log());
      }
      else {
        custom_generators.insert(copy_apply(op, val.cluster));
      }
    }
    return custom_generators;
  }

  // --- template<typename PhenomenalType> ClusterEquivalenceParser ---

  template<typename PhenomenalType>
  ClusterEquivalenceParser<PhenomenalType>::ClusterEquivalenceParser(
    const Supercell &_scel,
    jsonParser &_input,
    fs::path _path,
    bool _required):
    KwargsParser(_input, _path, _required),
    primclex(_scel.primclex()),
    scel(_scel) {

    // try to read phenom
    try {
      phenom = notstd::make_unique<PhenomenalType>(
                 jsonConstructor<PhenomenalType>::from_json(self, primclex));
    }
    catch(std::exception &e) {
      error.insert(std::string("Error: ") + "Could not read phenomal cluster: '" + e.what() + "'");
    }

    // try to read equivalence type, run checks, and generate necessary data
    self.get_else(equiv_type, traits<EQUIVALENCE_TYPE>::name, EQUIVALENCE_TYPE::PRIM);
    if(equiv_type == EQUIVALENCE_TYPE::PRIM) {
      _init_prim_equivalence();
    }
    else if(equiv_type == EQUIVALENCE_TYPE::SCEL) {
      _init_scel_equivalence();
    }
    else if(equiv_type == EQUIVALENCE_TYPE::CONFIG) {
      _init_config_equivalence();
    }
  }

  /// \brief Check if test is equivalent to phenom
  ///
  /// \returns (phenom.apply_sym(op) is_equivalent to test, op)
  ///
  /// For equiv_type == EQUIVALENCE_TYPE::PRIM, use prim factor group
  /// For equiv_type == EQUIVALENCE_TYPE::SCEL, use factor group of Supercell
  /// For equiv_type == EQUIVALENCE_TYPE::CONFIG, use factor group of Configuration named
  ///   by 'configname' (must have same Supercell as given to constructor)
  ///
  template<typename PhenomenalType>
  std::pair<bool, SymOp> ClusterEquivalenceParser<PhenomenalType>::is_equivalent(
    const PhenomenalType &test) const {
    if(equiv_type == EQUIVALENCE_TYPE::PRIM) {
      return _is_equivalent(test, primclex.prim().factor_group(), *prim_sym_compare);
    }
    else if(equiv_type == EQUIVALENCE_TYPE::SCEL) {
      // use prim_sym_compare because any translation is OK
      return _is_equivalent(test, scel.factor_group(), *prim_sym_compare);
    }
    else if(equiv_type == EQUIVALENCE_TYPE::CONFIG) {
      return _is_equivalent(test, *config_sym_compare, config_fg.begin(), config_fg.end());
    }
    throw std::runtime_error("Error: Unknown error in ClusterEquivalenceParser::is_equivalent");
  }

  template<typename PhenomenalType>
  void ClusterEquivalenceParser<PhenomenalType>::_init_prim_equivalence() {
    prim_sym_compare = notstd::make_unique<PrimPeriodicSymCompare<PhenomenalType>>(primclex);
    if(phenom) {
      *phenom = prim_sym_compare->prepare(*phenom);
    }
    if(self.find("configname") != self.end()) {
      warning.insert("Warning: Ignoring unnecessary 'configname' option");
    }
  }

  template<typename PhenomenalType>
  void ClusterEquivalenceParser<PhenomenalType>::_init_scel_equivalence() {
    prim_sym_compare = notstd::make_unique<PrimPeriodicSymCompare<PhenomenalType>>(primclex);
    scel_sym_compare = notstd::make_unique<ScelPeriodicSymCompare<PhenomenalType>>(scel);
    if(phenom) {
      *phenom = scel_sym_compare->prepare(*phenom);
    }
    if(self.find("configname") != self.end()) {
      warning.insert("Warning: Ignoring unnecessary 'configname' option");
    }
  }

  template<typename PhenomenalType>
  void ClusterEquivalenceParser<PhenomenalType>::_init_config_equivalence() {
    auto configname_ptr = require<std::string>("configname");
    if(configname_ptr) {
      configname = *configname_ptr;
      if(Configuration::split_name(configname).first != scel.name()) {
        error.insert(std::string("Error:") + " Checking custom local clusters. "
                     + "The configuration specified for equivalence checking, '"
                     + configname + "' does that have the correct supercell, '" + scel.name() + "'");
      }
      else {
        auto config_it = primclex.db<Configuration>().find(configname);
        config_fg = config_it->factor_group();
        config_sym_compare = notstd::make_unique<ScelPeriodicSymCompare<PhenomenalType>>(
                               config_it->supercell());

        if(phenom) {
          *phenom = config_sym_compare->prepare(*phenom);
        }
      }
    }
    if(self.find("scelname") != self.end()) {
      warning.insert("Warning: Ignoring unnecessary 'scelname' option");
    }
  }

  template<typename PhenomenalType>
  template<typename SymCompareType>
  std::pair<bool, SymOp> ClusterEquivalenceParser<PhenomenalType>::_is_equivalent(
    const PhenomenalType &test,
    const SymGroup &g,
    const SymCompareType &sym_compare) const {
    PhenomenalType _test = sym_compare.prepare(test);
    for(const auto &op : g) {
      auto tmp = sym_compare.prepare(copy_apply(op, *phenom));
      if(sym_compare.equal(tmp, _test)) {
        return std::make_pair(true, sym_compare.spatial_transform() * op);
      }
    }
    return std::make_pair(false, g[0]);
  }

  template<typename PhenomenalType>
  template<typename SymCompareType>
  std::pair<bool, SymOp> ClusterEquivalenceParser<PhenomenalType>::_is_equivalent(
    const PhenomenalType &test,
    const SymCompareType &sym_compare,
    PermuteIterator begin,
    PermuteIterator end) const {
    PhenomenalType _test = sym_compare.prepare(test);
    for(auto it = begin; it != end; ++it) {
      auto tmp = sym_compare.prepare(copy_apply(it, *phenom));
      if(sym_compare.equal(tmp, _test)) {
        return std::make_pair(true, sym_compare.spatial_transform() * it.sym_op());
      }
    }
    return std::make_pair(false, begin.sym_op());
  }

  template<typename PhenomenalType>
  template<typename SymCompareType, typename PermuteIteratorIt>
  std::pair<bool, SymOp> ClusterEquivalenceParser<PhenomenalType>::_is_equivalent(
    const PhenomenalType &test,
    const SymCompareType &sym_compare,
    PermuteIteratorIt begin,
    PermuteIteratorIt end) const {
    PhenomenalType _test = sym_compare.prepare(test);
    for(auto it = begin; it != end; ++it) {
      auto tmp = sym_compare.prepare(copy_apply(*it, *phenom));
      if(sym_compare.equal(tmp, _test)) {
        return std::make_pair(true, sym_compare.spatial_transform() * it->sym_op());
      }
    }
    return std::make_pair(false, begin->sym_op());
  }


  // --- template<typename PhenomenalType> LocalClustersByMaxLength ---

  template<typename PhenomenalType>
  CustomLocalClustersByMaxLength<PhenomenalType>::CustomLocalClustersByMaxLength(
    const Supercell &_scel,
    jsonParser &_input,
    fs::path _path,
    bool _required):
    InputParser(_input, _path, _required) {

    if(exists()) {
      if(!self.is_array()) {
        error.insert(std::string("Error: '") + name() + "' is not a JSON array");
        return;
      }

      for(Index i = 0; i < self.size(); ++i) {
        data.push_back(CustomSpecs());

        fs::path p_base = boost::lexical_cast<std::string>(i);

        fs::path p_phenom = p_base / "phenomenal";
        data.back().phenom = std::make_shared<ClusterEquivalenceParser<PhenomenalType>>(
                               _scel, input, relpath(p_phenom), true);
        this->kwargs[p_phenom] = data.back().phenom;

        fs::path p_branch_specs = p_base / "orbit_branch_specs";
        data.back().orbit_branch_specs = std::make_shared<LocalOrbitBranchSpecsParser>(
                                           input, relpath(p_branch_specs), false);
        this->kwargs[p_branch_specs] = data.back().orbit_branch_specs;

        fs::path p_orbit_specs = p_base / "orbit_specs";
        data.back().orbit_specs = std::make_shared<LocalOrbitSpecsParser>(
                                    _scel.primclex(), input, relpath(p_orbit_specs), false);
        this->kwargs[p_orbit_specs] = data.back().orbit_specs;
      }
    }
  }

  /// Find if phenom is equivalent to one of the custom phenomenal clusters
  ///
  /// \returns pair of: iterator to CustomSpecs for equivalent phenomenal cluster (if found),
  /// and the op such that custom_phenom.apply_sym(op) == phenom
  ///
  /// - If result.first == custom.end(), use standard specs; else use custom specs
  ///
  template<typename PhenomenalType>
  std::pair<typename CustomLocalClustersByMaxLength<PhenomenalType>::custom_specs_iterator, SymOp>
  CustomLocalClustersByMaxLength<PhenomenalType>::find(const PhenomenalType &phenom) const {
    for(auto it = data.begin(); it != data.end(); ++it) {
      auto res = it->phenom->is_equivalent(phenom);
      if(res.first) {
        return std::make_pair(it, res.second);
      }
    }
    return std::make_pair(data.end(), SymOp());
  }


  // --- template<typename PhenomenalType> LocalClustersByMaxLength ---

  template<typename PhenomenalType>
  LocalClustersByMaxLength<PhenomenalType>::LocalClustersByMaxLength(
    const Supercell &_scel,
    jsonParser &_input,
    fs::path _path,
    bool _required):
    InputParser(_input, _path, _required) {

    // Add parsers for standard clusters
    fs::path p = fs::path("standard") / "orbit_branch_specs";
    standard = std::make_shared<LocalOrbitBranchSpecsParser>(input, relpath(p), false);
    this->kwargs[p] = standard;

    // Add parsers for custom clusters
    p = fs::path("custom");
    custom = std::make_shared<CustomLocalClustersByMaxLength<PhenomenalType>>(
               _scel, input, relpath(p), false);
    this->kwargs[p] = custom;
    if(!(standard->exists() || custom->exists()))
      error.insert("Error: Either 'standard' : {'orbit_branch_specs' : {<...>}} must exists or 'custom' : {<...>}.");
  }

  /// Find if phenom is equivalent to one of the custom phenomenal clusters
  ///
  /// \returns pair of: iterator to CustomSpecs for equivalent phenomenal cluster (if found),
  /// and the op such that custom_phenom.apply_sym(op) == phenom
  ///
  /// - If result.first == custom.end(), use standard specs; else use custom specs
  ///
  template<typename PhenomenalType>
  std::pair<typename LocalClustersByMaxLength<PhenomenalType>::custom_specs_iterator, SymOp>
  LocalClustersByMaxLength<PhenomenalType>::find(const PhenomenalType &phenom) const {
    return custom->find(phenom);
  }

  template<typename PhenomenalType>
  bool LocalClustersByMaxLength<PhenomenalType>::max_length_including_phenomenal(custom_specs_iterator it) const {
    if(it == custom->data.end()) {
      return standard->max_length_including_phenomenal();
    }
    return it->orbit_branch_specs->max_length_including_phenomenal();
  }

  template<typename PhenomenalType>
  int LocalClustersByMaxLength<PhenomenalType>::max_branch(custom_specs_iterator it) const {
    if(it == custom->data.end()) {
      return standard->max_branch;
    }
    return it->orbit_branch_specs->max_branch;
  }

  template<typename PhenomenalType>
  double LocalClustersByMaxLength<PhenomenalType>::max_length(custom_specs_iterator it, int branch_i) const {
    if(it == custom->data.end()) {
      return standard->max_length(branch_i);
    }
    return it->orbit_branch_specs->max_length(branch_i);
  }

  template<typename PhenomenalType>
  double LocalClustersByMaxLength<PhenomenalType>::cutoff_radius(custom_specs_iterator it, int branch_i) const {
    if(it == custom->data.end()) {
      return standard->cutoff_radius(branch_i);
    }
    return it->orbit_branch_specs->cutoff_radius(branch_i);
  }

  // *INDENT-OFF*

  /// Get custom local cluster generators
  ///
  /// \throws std::invalid_argument if find_res.first == custom.end()
  template<typename PhenomenalType>
  template<typename OrbitType>
  OrbitGenerators<OrbitType>&
  LocalClustersByMaxLength<PhenomenalType>::insert_custom_generators(
      std::pair<custom_specs_iterator, SymOp> find_res,
      OrbitGenerators<OrbitType>& custom_generators) const {
    if(find_res.first == custom->data.end()) {
      throw std::invalid_argument("Error: No custom generators for standard local orbit specs");
    }
    return find_res.first->orbit_specs->insert_custom_generators(find_res.second, custom_generators);
  }

  // *INDENT-ON*
}

#endif
