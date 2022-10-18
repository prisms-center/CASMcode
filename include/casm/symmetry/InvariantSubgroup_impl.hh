#ifndef CASM_InvariantSubgroup_impl
#define CASM_InvariantSubgroup_impl

#include "casm/symmetry/InvariantSubgroup.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/symmetry/SymGroup.hh"

namespace CASM {

// using xtal::Coordinate;

/// Construct the subgroup that leaves an element unchanged
///
/// Includes translations as determined from 'sym_compare'
///
/// Implementation is equivalent to:
/// \code
/// Element e {sym_compare.prepare(element)};
/// for(auto group_it=group_begin; group_it!=group_end; ++group_it) {
///   if(sym_compare.equal(e, sym_compare.prepare(sym_compare.copy_apply(*it,
///   e)))) {
///     *result++ = sym_compare.spatial_transform() * (*it);
///   }
/// }
/// return result;
/// \endcode
template <typename Element, typename OpIterator, typename SymCompareType,
          typename OpOutputIterator>
OpOutputIterator make_invariant_subgroup(Element const &element,
                                         OpIterator group_begin,
                                         OpIterator group_end,
                                         SymCompareType const &sym_compare,
                                         OpOutputIterator result) {
  Element e{sym_compare.prepare(element)};
  for (auto group_it = group_begin; group_it != group_end; ++group_it) {
    if (sym_compare.equal(
            e, sym_compare.prepare(sym_compare.copy_apply(*group_it, e)))) {
      *result++ = sym_compare.spatial_transform() * (*group_it);
    }
  }
  return result;
}

/// Construct the subgroup (PermuteIterator) that leaves an element unchanged
///
/// Implementation is equivalent to:
/// \code
/// Element e {sym_compare.prepare(element)};
/// for(auto permute_it=permute_begin; permute_it!=permute_end; ++permute_it) {
///   if(sym_compare.equal(e,
///   sym_compare.prepare(sym_compare.copy_apply(permute_it->sym_op(), e)))) {
///     *result++ = permute_it;
///   }
/// }
/// return result;
/// \endcode
template <typename Element, typename SymCompareType, typename OutputIterator>
OutputIterator make_invariant_subgroup(Element const &element,
                                       PermuteIterator permute_begin,
                                       PermuteIterator permute_end,
                                       SymCompareType const &sym_compare,
                                       OutputIterator result) {
  Element e{sym_compare.prepare(element)};
  for (auto permute_it = permute_begin; permute_it != permute_end;
       ++permute_it) {
    if (sym_compare.equal(e, sym_compare.prepare(sym_compare.copy_apply(
                                 permute_it->sym_op(), e)))) {
      *result++ = permute_it;
    }
  }
  return result;
}

/// \brief Construct the subgroup that leaves an element unchanged
///
/// Includes translations as determined from 'sym_compare'
///
/// Implementation is equivalent to:
/// \code
/// Element e {sym_compare.prepare(element)};
/// SymGroup result = generating_grp;
/// result.clear();
/// for(const auto &op : generating_grp) {
///   if(sym_compare.equal(e, sym_compare.prepare(sym_compare.copy_apply(op,
///   e)))) {
///    result.push_back(sym_compare.spatial_transform()*op);
///   }
/// }
/// return result;
/// \endcode
template <typename Element, typename SymCompareType>
SymGroup make_invariant_subgroup(const Element &element,
                                 const SymGroup &generating_grp,
                                 const SymCompareType &sym_compare) {
  if (element.size() == 0) {
    return generating_grp;
  }
  SymGroup result = generating_grp;
  result.clear();
  make_invariant_subgroup(element, generating_grp.begin(), generating_grp.end(),
                          sym_compare, std::back_inserter(result));
  if (result[0].index() != 0) {
    throw std::runtime_error(
        "Error in make_invariant_subgroup (0): First element is not identity.");
  }
  return result;
}

/// \brief Construct the subgroup that leaves an element of the orbit unchanged
///
/// Does not include translations as determined from 'sym_compare'
///
/// The equivalence map is:
///   element(i) compares equivalent to
///   prototype().copy_apply(equivalence_map[i][j]) for all j
///
/// The subgroup that leaves element 'i' of the orbit unchanged is, for all j:
///   eq_map[i][0]*eq_map[0][j]*eq_map[i][0].inverse()
///
template <typename OrbitType>
SymGroup make_invariant_subgroup(const OrbitType &orbit, Index element_index) {
  SymGroup result;
  const auto &map = orbit.equivalence_map();
  result.set_lattice(map[0][0].master_group().lattice());
  for (Index i = 0; i < orbit.equivalence_map()[0].size(); ++i) {
    result.push_back(map[element_index][0] * map[0][i] *
                     map[element_index][0].inverse());
  }
  if (result[0].index() != 0) {
    throw std::runtime_error(
        "Error in make_invariant_subgroup (1): First element is not identity.");
  }
  return result;
}

}  // namespace CASM

#endif
