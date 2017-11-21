#ifndef CASM_ClusterSpecsParser_impl
#define CASM_ClusterSpecsParser_impl

#include "casm/clusterography/ClusterSpecsParser.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/CASM_global_definitions.hh"
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
    if(prev_it == self_it->end()) {
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

  // --- template<typename PhenomenalType> LocalOrbitSpecsParser ---

  template<typename PhenomenalType>
  LocalOrbitSpecsParser<PhenomenalType>::LocalOrbitSpecsParser(
    jsonParser &_input,
    fs::path _path,
    bool _required) {
    //...todo...
  }

  /// For all custom clusters, insert 'op*prototype' into custom_generators
  ///
  /// \note Use ClusterEquivalenceParser to determine if custom clusters apply
  /// to a given cluster, and determine 'op'.
  ///
  template<typename PhenomenalType>
  const OrbitGenerators<LocalOrbitSpecsParser<PhenomenalType>::OrbitType> &
  LocalOrbitSpecsParser<PhenomenalType>::custom_generators(
    const SymOp &op,
    const SymGroup &generating_grp,
    const OrbitType::SymCompareType &sym_compare) const {
    // ...todo...
  }

  // --- template<typename PhenomenalType> ClusterEquivalenceParser ---

  template<typename PhenomenalType>
  ClusterEquivalenceParser<PhenomenalType>::ClusterEquivalenceParser(
    const PrimClex &_primclex,
    jsonParser &_input,
    fs::path _path,
    bool _required):
    primclex(_primclex),
    scel(nullptr) {

    // try to read phenom
    try {
      phenom = notstd::make_unique<PhenomenalType>(
                 jsonConstructor<PhenomenalType>::from_json(*self_it, primclex));
    }
    catch(std::exception &e) {
      error.insert("Error: Could not read phenomal cluster");
    }

    // try to read equivalence type, run checks, and generate necessary data
    self_it->get_else(equiv_type, traits<EQUIVALENCE_TYPE>::name, EQUIVALENCE_TYPE::PRIM);
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
  /// For equiv_type == EQUIVALENCE_TYPE::SCEL, use factor group of Supercell named by 'scelname'
  /// For equiv_type == EQUIVALENCE_TYPE::CONFIG, use factor group of Configuration named by 'configname'
  ///
  template<typename PhenomenalType>
  std::pair<bool, SymOp> ClusterEquivalenceParser<PhenomenalType>::is_equivalent(const PhenomenalType &test) const {
    if(equiv_type == EQUIVALENCE_TYPE::PRIM) {
      return _is_equivalent(test, primclex.prim().factor_group(), *prim_sym_compare);
    }
    else if(equiv_type == EQUIVALENCE_TYPE::SCEL) {
      return _is_equivalent(test, *scel_sym_compare, scel->permute_begin(), scel->permute_end());
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
      *phenom = scel_sym_compare.prepare(*phenom);
    }
    if(self_it->find("scelname")) {
      warning.insert("Warning: Ignoring unnecessary 'scelname' option");
    }
    if(self_it->find("configname")) {
      warning.insert("Warning: Ignoring unnecessary 'configname' option");
    }
  }

  template<typename PhenomenalType>
  void ClusterEquivalenceParser<PhenomenalType>::_init_scel_equivalence() {
    auto scelname = require<std::string>("scelname");
    if(scelname) {
      auto scel_it = primclex.db<Supercell>().find(*scelname);
      scel = &(*scel_it);
      scel_sym_compare = notstd::make_unique<ScelPeriodicSymCompare<PhenomenalType>>(*scel);

      if(phenom) {
        *phenom = scel_sym_compare->prepare(*phenom);
      }
    }
    if(self_it->find("configname")) {
      warning.insert("Warning: Ignoring unnecessary 'configname' option");
    }
  }

  template<typename PhenomenalType>
  void ClusterEquivalenceParser<PhenomenalType>::_init_config_equivalence() {
    auto configname = require<std::string>("configname");
    if(configname) {
      auto scel_it = primclex.db<Configuration>().find(*configname);
      scel = &(*scel_it);
      config_sym_compare = notstd::make_unique<ScelPeriodicSymCompare<PhenomenalType>>(*scel);

      if(phenom) {
        *phenom = config_sym_compare->prepare(*phenom);
      }
    }
    if(self_it->find("scelname")) {
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
      auto tmp = sym_compare.prepare(copy_apply(op, phenom));
      if(sym_compare.equal(tmp, _test)) {
        return std::make_pair(true, op);
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
      auto tmp = sym_compare.prepare(copy_apply(it, phenom));
      if(sym_compare.equal(tmp, _test)) {
        return std::make_pair(true, it.sym_op());
      }
    }
    return std::make_pair(false, begin);
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
      auto tmp = sym_compare.prepare(copy_apply(*it, phenom));
      if(sym_compare.equal(tmp, _test)) {
        return std::make_pair(true, it->sym_op());
      }
    }
    return std::make_pair(false, begin->sym_op());
  }


  // --- template<typename PhenomenalType> LocalClustersByMaxLength ---

  template<typename PhenomenalType>
  LocalClustersByMaxLength<PhenomenalType>::LocalClustersByMaxLength(
    const PrimClex &_primclex,
    jsonParser &_input,
    fs::path _path,
    bool _required) {
    // ...todo...
  }

  /// Find if phenom is equivalent to one of the custom phenomenal clusters
  ///
  /// \returns pair of: iterator to CustomSpecs for equivalent phenomenal cluster (if found),
  /// and the op such that custom_phenom.apply_sym(op) == phenom
  ///
  /// - If result.first == custom.end(), use standard specs; else use custom specs
  ///
  template<typename PhenomenalType>
  std::pair<custom_specs_iterator, SymOp>
  LocalClustersByMaxLength<PhenomenalType>::find(const PhenomenalType &phenom) const {
    for(const auto &custom_specs : custom) {
      auto res = custom_specs.phenom.is_equivalent(phenom);
      if(res.first) {
        return res;
      }
    }
    return std::make_pair(custom.end(), SymOp());
  }

  template<typename PhenomenalType>
  bool LocalClustersByMaxLength<PhenomenalType>::max_length_including_phenomenal(custom_specs_iterator it) const {
    if(it == custom.end()) {
      return standard.max_length_including_phenomenal();
    }
    return it->orbit_branch_specs.max_length_including_phenomenal();
  }

  template<typename PhenomenalType>
  int LocalClustersByMaxLength<PhenomenalType>::max_branch(custom_specs_iterator it) const {
    if(it == custom.end()) {
      return standard.max_branch();
    }
    return it->orbit_branch_specs.max_branch();
  }

  template<typename PhenomenalType>
  double LocalClustersByMaxLength<PhenomenalType>::max_length(custom_specs_iterator it, int branch_i) const {
    if(it == custom.end()) {
      return standard.max_length(branch_i);
    }
    return it->orbit_branch_specs.max_length(branch_i);
  }

  template<typename PhenomenalType>
  double LocalClustersByMaxLength<PhenomenalType>::cutoff_radius(custom_specs_iterator it, int branch_i) const {
    if(it == custom.end()) {
      return standard.cutoff_radius(branch_i);
    }
    return it->orbit_branch_specs.cutoff_radius(branch_i);
  }

  // *INDENT-OFF*

  /// Get custom local cluster generators
  ///
  /// \throws std::invalid_argument if find_res.first == custom.end()
  template<typename PhenomenalType>
  template<typename SymCompareType>
  const OrbitGenerators<Orbit<IntegralCluster, SymCompareType>> &
  LocalClustersByMaxLength<PhenomenalType>::custom_generators(
      std::pair<custom_specs_iterator, SymOp> find_res,
      const SymGroup &generating_grp,
      const SymCompareType &sym_compare) const {
    if(find_res.first == custom.end()) {
      throw std::invalid_argument("Error: No custom generators for standard local orbit specs");
    }
    return find_res.first->orbit_specs.custom_generators(find_res.second, generating_grp, sym_compare);
  }

  // *INDENT-ON*
}

#endif
