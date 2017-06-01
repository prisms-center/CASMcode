#ifndef CASM_InvariantSubgroup
#define CASM_InvariantSubgroup

#include <vector>
#include "casm/CASM_global_definitions.hh"

namespace CASM {

  class SymGroup;
  class PermuteIterator;
  class Supercell;

  /// \brief Construct the subgroup that leaves an element unchanged
  template<typename Element, typename SymCompareType>
  SymGroup make_invariant_subgroup(const Element &element,
                                   const SymGroup &generating_grp,
                                   const SymCompareType &sym_compare);

  /// \brief Construct the subgroup that leaves an element of the orbit unchanged
  template<typename OrbitType>
  SymGroup make_invariant_subgroup(const OrbitType &orbit, Index element_index = 0);

  /// \brief Construct the subgroup of permutations that leaves an element unchanged
  template<typename Element>
  std::vector<PermuteIterator> make_invariant_subgroup(
    const Element &element,
    const Supercell &scel);

}

#endif
