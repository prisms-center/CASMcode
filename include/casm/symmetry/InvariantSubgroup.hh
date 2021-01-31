#ifndef CASM_InvariantSubgroup
#define CASM_InvariantSubgroup

#include <vector>

#include "casm/global/definitions.hh"

namespace CASM {

class PermuteIterator;
class SymGroup;

/// Construct the subgroup that leaves an element unchanged
template <typename Element, typename OpIterator, typename SymCompareType,
          typename OpOutputIterator>
OpOutputIterator make_invariant_subgroup(Element const &element,
                                         OpIterator group_begin,
                                         OpIterator group_end,
                                         SymCompareType const &sym_compare,
                                         OpOutputIterator result);

/// Construct the subgroup (PermuteIterator) that leaves an element unchanged
template <typename Element, typename SymCompareType, typename OutputIterator>
OutputIterator make_invariant_subgroup(Element const &element,
                                       PermuteIterator permute_begin,
                                       PermuteIterator permute_end,
                                       SymCompareType const &sym_compare,
                                       OutputIterator result);

/// Construct the subgroup that leaves an element unchanged
template <typename Element, typename SymCompareType>
SymGroup make_invariant_subgroup(const Element &element,
                                 const SymGroup &generating_grp,
                                 const SymCompareType &sym_compare);

/// Construct the subgroup that leaves an element of the orbit unchanged
template <typename OrbitType>
SymGroup make_invariant_subgroup(const OrbitType &orbit,
                                 Index element_index = 0);

}  // namespace CASM

#endif
