#ifndef CASM_symmetry_SimpleOrbit_impl
#define CASM_symmetry_SimpleOrbit_impl

#include "casm/symmetry/SimpleOrbit.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {

  // --- template<typename _SymCompareType> class SimpleOrbitElementCompare ---

  template<typename _SymCompareType>
  SimpleOrbitElementCompare<_SymCompareType>::SimpleOrbitElementCompare(SymCompareType const &_sym_compare):
    sym_compare(_sym_compare) {}

  template<typename _SymCompareType>
  bool SimpleOrbitElementCompare<_SymCompareType>::operator()(const Element &A, const Element &B) const {
    return sym_compare.compare(A, B);
  }

  /// Construct a SimpleOrbit using a range of SymOp or PermuteIterator
  template<typename SymCompareType>
  template<typename GroupIterator>
  SimpleOrbit<SymCompareType>::SimpleOrbit(
    typename SymCompareType::Element const &_generating_element,
    GroupIterator _group_begin,
    GroupIterator _group_end,
    SymCompareType const &_sym_compare):
    m_sym_compare(_sym_compare),
    m_invariants(_sym_compare.make_invariants(_generating_element)) {

    std::set<Element, SimpleOrbitElementCompare<SymCompareType> > tmp {
      SimpleOrbitElementCompare<SymCompareType>{m_sym_compare}};
    adapter::Adapter<SymOp, typename GroupIterator::value_type> to_symop;
    for(auto it = _group_begin; it != _group_begin; ++it) {
      tmp.insert(m_sym_compare.prepare(m_sym_compare.copy_apply(to_symop(*it), _generating_element)));
    }

    std::copy(tmp.begin(), tmp.end(), std::back_inserter(m_element));
  }

  /// Compare orbits, using SymCompareType::inter_orbit_compare
  template<typename SymCompareType>
  bool SimpleOrbit<SymCompareType>::operator<(const SimpleOrbit &B) const {
    return m_sym_compare.inter_orbit_compare(prototype(), invariants(), B.prototype(), B.invariants());
  }

}

#endif
