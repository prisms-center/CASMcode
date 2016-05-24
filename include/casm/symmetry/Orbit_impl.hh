#ifndef CASM_Orbit_impl
#define CASM_Orbit_impl

#include "casm/symmetry/Orbit.hh"

namespace CASM {

  /* -- Orbit Definitions ------------------------------------- */

  /// \brief Construct an Orbit from a generating_element Element, using provided Group
  ///
  /// \param generating_element The element used to generate equivalents
  /// \param generating_group The group used for generating the orbit
  /// \param sym_compare Binary functor that implements symmetry properties
  ///
  template<typename Element>
  Orbit<Element>::Orbit(Element generating_element,
                        const SymGroup &generating_group,
                        const SymCompare<Element> &sym_compare) :
    m_sym_compare(sym_compare) {
    _construct(generating_element, generating_group.begin(), generating_group.end());
  }

  /// \brief Construct an Orbit from a generating_element Element, using provided symmetry rep
  ///
  /// \param generating_element The element used to generate equivalents
  /// \param begin, end Range of SymOp applied to Element to generate the orbit
  /// \param sym_compare Binary functor that implements symmetry properties
  ///
  /// - iterators must be multi-pass
  template<typename Element>
  template<typename SympOpIterator>
  Orbit<Element>::Orbit(Element generating_element,
                        SympOpIterator begin,
                        SympOpIterator end,
                        const SymCompare<Element> &sym_compare) :
    m_sym_compare(sym_compare) {
    _construct(generating_element, begin, end);
  }

  /// \brief Construct an Orbit from a generating_element Element, using provided symmetry rep
  ///
  /// \param generating_element The element used to generate equivalents
  /// \param begin, end Range of SymOp applied to Element to generate the orbit
  ///
  /// - iterators must be multi-pass
  template<typename Element>
  template<typename SympOpIterator>
  void Orbit<Element>::_construct(Element generating_element,
                                  SympOpIterator begin,
                                  SympOpIterator end) {

    // add all op*generating_element
    for(auto it = begin; it != end; ++it) {
      m_element.push_back(m_sym_compare->prepare(copy_apply(*it, generating_element)));
    }

    // sort
    std::sort(m_element.begin(), m_element.end(), *m_sym_compare);

    // keep uniques
    auto last = std::unique(
                  m_element.begin(),
                  m_element.end(),
    [](const Element & A, const Element & B) {
      return m_sym_compare->intra_orbit_equal(A, B);
    });
    m_element.erase(last, m_element.end());

    // create equivalence map
    for(auto it = begin; it != end; ++it) {
      auto equiv = m_sym_compare->prepare(copy_apply(*it, prototype()));
      equivalence_map[find_index(m_element, equiv)].push_back(*it);
    }
  }

  /// \brief Apply symmetry to Orbit
  template<typename Element>
  Orbit<Element> &Orbit<Element>::apply_sym(const SymOp &op) {

    // transform elements
    for(auto it = m_element.begin(); it != m_element.end(); ++it) {
      it->apply_sym(op);
    }

    // transform equivalence map: std::vector<std::vector<SymOp> >
    for(auto it = m_equivalence_map.begin(); it != m_equivalence_map.end(); ++it) {
      for(auto op_it = it->begin(); op_it != it->end(); ++op_it) {
        op_it->apply_sym(op);
      }
    }

    // transform sym_compare functor
    m_sym_compare->apply_sym(op);

    return *this;
  }


}

#endif