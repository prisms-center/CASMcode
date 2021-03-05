#ifndef CASM_symmetry_SimpleOrbit_impl
#define CASM_symmetry_SimpleOrbit_impl

#include <map>

#include "casm/global/errors.hh"
#include "casm/symmetry/SimpleOrbit.hh"

namespace CASM {

/// Construct a SimpleOrbit using a range of SymOp or PermuteIterator
template <typename SymCompareType>
template <typename GroupIterator>
SimpleOrbit<SymCompareType>::SimpleOrbit(
    typename SymCompareType::Element const &_generating_element,
    GroupIterator _group_begin, GroupIterator _group_end,
    SymCompareType const &_sym_compare)
    : m_sym_compare(_sym_compare),
      m_invariants(_sym_compare.make_invariants(_generating_element)) {
  auto compare = [&](const Element &A, const Element &B) {
    return m_sym_compare.compare(A, B);
  };

  typedef typename GroupIterator::value_type SymOpRepType;
  std::map<Element, std::set<SymOpRepType>, decltype(compare)> tmp{compare};
  for (auto it = _group_begin; it != _group_end; ++it) {
    tmp[m_sym_compare.prepare(
            m_sym_compare.copy_apply(*it, _generating_element))]
        .insert(*it);
  }

  // sanity check equivalence map is rectangular
  for (auto const &pair : tmp) {
    if (tmp.begin()->second.size() != pair.second.size()) {
      throw libcasm_runtime_error(
          "Error in SimpleOrbit constructor: equivalence map is not "
          "rectangular");
    }
  }

  m_element.reserve(tmp.size());
  for (auto const &pair : tmp) {
    m_element.push_back(pair.first);
  }
}

/// Compare orbits, using SymCompareType::inter_orbit_compare
template <typename SymCompareType>
bool SimpleOrbit<SymCompareType>::operator<(const SimpleOrbit &B) const {
  return m_sym_compare.inter_orbit_compare(prototype(), invariants(),
                                           B.prototype(), B.invariants());
}

}  // namespace CASM

#endif
