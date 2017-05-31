#ifndef CASM_SubOrbits_impl
#define CASM_SubOrbits_impl

#include "casm/symmetry/InvariantSubgroup_impl.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/symmetry/SubOrbits.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to group -> subgroup symmetry breaking
  template<typename Element, typename ElementOutputIterator>
  ElementOutputIterator make_suborbit_generators(
    const Element &element,
    const SymGroup &invariant_subgroup,
    const SymGroup &group,
    const SymGroup &subgroup,
    ElementOutputIterator result) {

    // find "max" SymOp in each coset of the subgroup,
    //   excluding those cosets that cause duplicate generating elements
    for(Index i = 0; i < group.size(); ++i) {
      const SymOp &test_op = group[i];
      auto lambda = [&](const SymOp & op) {
        for(const auto &el_op : invariant_subgroup) {
          if(test_op.index() < (op * test_op * el_op).index()) {
            return true;
          }
        }
        return false;
      };
      // if test_op is max
      if(std::none_of(subgroup.begin(), subgroup.end(), lambda)) {
        // apply to element and construct suborbit generator
        *result++ = copy_apply(test_op, element);
      }
    }
    return result;
  }

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to group -> subgroup symmetry breaking
  template<typename Element, typename SymCompareType, typename ElementOutputIterator>
  ElementOutputIterator make_suborbit_generators(
    const Element &element,
    const SymCompareType &sym_compare,
    const SymGroup &group,
    const SymGroup &subgroup,
    ElementOutputIterator result) {

    SymGroup invariant_subgroup = make_invariant_subgroup(element, group, sym_compare);
    return make_suborbit_generators(element, invariant_subgroup, group, subgroup, result);
  }

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to group -> subgroup symmetry breaking
  template<typename OrbitType, typename ElementOutputIterator>
  ElementOutputIterator make_suborbit_generators(
    const OrbitType &orbit,
    const SymGroup &group,
    const SymGroup &subgroup,
    ElementOutputIterator result) {

    SymGroup invariant_subgroup = make_invariant_subgroup(orbit);
    return make_suborbit_generators(orbit.prototype(), invariant_subgroup, group, subgroup, result);
  }

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to group -> subgroup symmetry breaking
  template<typename Element, typename ElementOutputIterator, typename PermuteIteratorIt>
  ElementOutputIterator make_suborbit_generators(
    const Element &element,
    const Supercell &scel,
    PermuteIteratorIt subgroup_begin,
    PermuteIteratorIt subgroup_end,
    ElementOutputIterator result) {

    const PrimGrid &prim_grid = scel.prim_grid();

    std::vector<PermuteIterator> scel_inv_group = make_invariant_subgroup(element, scel);

    // find "max" permute in each coset of the subgroup [begin, end),
    //   excluding those cosets that cause duplicate generating elements
    auto test_it = scel.permute_begin();
    auto test_end = scel.permute_end();
    for(; test_it != test_end; ++test_it) {
      auto lambda = [&](const PermuteIterator & permute_it) {
        for(auto it = scel_inv_group.begin(); it != scel_inv_group.end(); ++it) {
          if(test_it < permute_it * (test_it * (*it))) {
            return true;
          }
        }
        return false;
      };

      // if test_it is max
      if(std::none_of(subgroup_begin, subgroup_end, lambda)) {
        // apply to prototype and construct suborbit
        *result++ = copy_apply(test_it.sym_op(), element);
      }
    }

    return result;
  }

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to Prim Structure -> Configuration symmetry breaking
  ///
  /// Uses config.supercell() PermuteIterators, so slow if config is not primitive
  template<typename OrbitIterator, typename ElementOutputIterator>
  ElementOutputIterator make_suborbit_generators_slow(
    OrbitIterator begin,
    OrbitIterator end,
    const Configuration &config,
    ElementOutputIterator result) {

    typedef typename std::iterator_traits<OrbitIterator>::value_type OrbitType;
    typedef typename OrbitType::Element Element;

    const Structure &prim = config.prim();
    const Supercell &scel = config.supercell();
    std::vector<PermuteIterator> config_fg = config.factor_group();

    for(auto it = begin; it != end; ++it) {

      // get generating elements for prim->supercell orbit splitting
      std::vector<Element> scel_suborbit_generators;

      make_suborbit_generators(
        *it,
        prim.factor_group(),
        scel.factor_group(),
        std::back_inserter(scel_suborbit_generators));

      // get generating elements for supercell->config orbit splitting
      for(const auto &el : scel_suborbit_generators) {

        result = make_suborbit_generators(
                   el,
                   scel,
                   config_fg.begin(),
                   config_fg.end(),
                   result);
      }
    }
    return result;
  }

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to Prim Structure -> Configuration symmetry breaking
  ///
  /// Uses make_suborbit_generators_slow with primitive configuration, then
  /// splits again for case that the config supercell has lower symmetry than
  /// the primitive config supercell.
  template<typename OrbitIterator, typename ElementOutputIterator>
  ElementOutputIterator make_suborbit_generators(
    OrbitIterator begin,
    OrbitIterator end,
    const Configuration &config,
    ElementOutputIterator result) {

    typedef typename std::iterator_traits<OrbitIterator>::value_type OrbitType;
    typedef typename OrbitType::Element Element;

    Configuration prim_config = config.primitive().in_canonical_supercell();
    ScelPeriodicSymCompare<Element> sym_compare(
      config.supercell().prim_grid(),
      config.crystallography_tol());

    std::vector<Element> prim_config_suborbit_generators;

    make_suborbit_generators_slow(
      begin,
      end,
      prim_config,
      std::back_inserter(prim_config_suborbit_generators));

    for(const auto &el : prim_config_suborbit_generators) {
      result = make_suborbit_generators(
                 el,
                 sym_compare,
                 prim_config.supercell().factor_group(),
                 config.supercell().factor_group(),
                 result);
    }

    return result;
  }

}

#endif
