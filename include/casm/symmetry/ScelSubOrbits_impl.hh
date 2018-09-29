#ifndef CASM_ScelSubOrbits_impl
#define CASM_ScelSubOrbits_impl

#include "casm/symmetry/InvariantSubgroup_impl.hh"
#include "casm/symmetry/ScelSubOrbits.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/clex/Supercell.hh"

namespace CASM {

  // --- Prim -> Supercell ---

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to Prim Structure -> Supercell symmetry breaking
  template<typename Element, typename ElementOutputIterator, typename PermuteIteratorIt>
  ElementOutputIterator make_suborbit_generators(
    const Element &element,
    const Supercell &scel,
    PermuteIteratorIt subgroup_begin,
    PermuteIteratorIt subgroup_end,
    ElementOutputIterator result) {

    if(&scel.prim_grid() != &subgroup_begin->prim_grid()) {
      throw std::runtime_error("Error: PermuteIterator supercell mismatch.");
    }

    std::vector<PermuteIterator> scel_inv_group = make_invariant_subgroup(element, scel);

    // find "max" permute in each coset of the subgroup [begin, end),
    //   excluding those cosets that cause duplicate generating elements
    auto test_it = scel.sym_info().permute_begin();
    auto test_end = scel.sym_info().permute_end();
    for(; test_it != test_end; ++test_it) {
      auto lambda = [&](const PermuteIterator & permute_it) {
        for(const auto &el_it : scel_inv_group) {
          if(test_it < (permute_it * test_it * el_it)) {
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
  /// corresponding to Prim Structure -> Supercell symmetry breaking
  template<typename Element, typename ElementOutputIterator, typename PermuteIteratorIt>
  ElementOutputIterator make_suborbit_generators(
    const Element &element,
    const Supercell &scel,
    PermuteIteratorIt group_begin,
    PermuteIteratorIt group_end,
    PermuteIteratorIt subgroup_begin,
    PermuteIteratorIt subgroup_end,
    ElementOutputIterator result) {

    if(&group_begin->prim_grid() != &subgroup_begin->prim_grid()) {
      throw std::runtime_error("Error: PermuteIterator supercell mismatch.");
    }

    std::vector<PermuteIterator> scel_inv_group = make_invariant_subgroup(
                                                    element, scel, group_begin, group_end);

    // find "max" permute in each coset of the subgroup [begin, end),
    //   excluding those cosets that cause duplicate generating elements
    auto test_it = group_begin;
    auto test_end = group_end;
    for(; test_it != test_end; ++test_it) {
      auto lambda = [&](const PermuteIterator & permute_it) {
        for(const auto &el_it : scel_inv_group) {
          if((*test_it) < (permute_it * (*test_it) * el_it)) {
            return true;
          }
        }
        return false;
      };

      // if test_it is max
      if(std::none_of(subgroup_begin, subgroup_end, lambda)) {
        // apply to prototype and construct suborbit
        *result++ = copy_apply(test_it->sym_op(), element);
      }
    }

    return result;
  }
}

#endif
