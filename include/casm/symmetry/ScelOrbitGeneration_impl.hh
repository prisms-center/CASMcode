#ifndef CASM_ScelOrbitGeneration_impl
#define CASM_ScelOrbitGeneration_impl

#include "casm/symmetry/ScelOrbitGeneration.hh"
#include "casm/clex/Supercell.hh"

namespace CASM {

  template<typename _ElementType>
  ScelCanonicalGenerator<_ElementType>::ScelCanonicalGenerator(const Supercell &_scel) :
    m_scel(&_scel),
    m_sym_compare(_scel.prim_grid(), _scel.crystallography_tol()) {}

  template<typename _ElementType>
  const Supercell &ScelCanonicalGenerator<_ElementType>::supercell() const {
    return *m_scel;
  }

  template<typename _ElementType>
  const typename ScelCanonicalGenerator<_ElementType>::SymCompareType &
  ScelCanonicalGenerator<_ElementType>::sym_compare() const {
    return m_sym_compare;
  }

  /// \brief Applies symmetry to return an equivalent Element in a canonical form
  ///
  /// - Use [supercell().permute_begin(), supercell().permute_end()) to canonicalize
  template<typename _ElementType>
  typename ScelCanonicalGenerator<_ElementType>::Element
  ScelCanonicalGenerator<_ElementType>::operator()(const Element &e) const {
    Element result = m_sym_compare.prepare(e);
    auto it = supercell().permute_begin();
    auto end = supercell().permute_end();
    while(it != end) {
      auto test = m_sym_compare.prepare(copy_apply(it, e));
      if(m_sym_compare.compare(result, test)) {
        result = test;
        m_to_canonical = it;
      }
      ++it;
    }
    return result;
  }

  /// \brief Applies symmetry to return an equivalent Element in a canonical form
  ///
  /// - For use with a container of PermuteIterator
  /// - Use [begin, end) to canonicalize
  template<typename _ElementType>
  template<typename PermuteIteratorIt>
  typename ScelCanonicalGenerator<_ElementType>::Element
  ScelCanonicalGenerator<_ElementType>::operator()(const Element &e, PermuteIteratorIt begin, PermuteIteratorIt end) const {
    Element result = m_sym_compare.prepare(e);
    auto it = begin;
    while(it != end) {
      auto test = m_sym_compare.prepare(copy_apply(*it, e));
      if(m_sym_compare.compare(result, test)) {
        result = test;
        m_to_canonical = *it;
      }
      ++it;
    }
    return result;
  }

  /// \brief After using call operator, this can be checked
  template<typename _ElementType>
  PermuteIterator ScelCanonicalGenerator<_ElementType>::to_canonical() const {
    return m_to_canonical;
  }

  /// \brief After using call operator, this can be checked
  template<typename _ElementType>
  PermuteIterator ScelCanonicalGenerator<_ElementType>::from_canonical() const {
    return to_canonical().inverse();
  }

}

#endif
