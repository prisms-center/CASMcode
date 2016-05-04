#ifndef CASM_Orbit
#define CASM_Orbit

#include <iostream>
#include <iomanip>

#include "casm/symmetry/UnitCellSiteSymOp.hh"
#include "casm/symmetry/FactorGroup.hh"

namespace casm {
  
  /* -- Orbit Declarations ------------------------------------- */
  
  /// \brief An Orbit of Element
  ///
  /// Provides prototype Element, orbit of equivalent elements, and equivalence 
  /// describing how the symmetry operations that map the prototype onto the
  /// equivalents.
  ///
  /// \ingroup Clusterography
  ///
  template<typename Element>
  class Orbit {
    
    public:
    
    typedef unsigned int size_type;
    typedef typename std::vector<Element>::value_type value_type;
    typedef typename std::vector<Element>::const_iterator const_iterator;
    
    Orbit() {}
    
    /// \brief Construct an Orbit from a generating_element Element, using provided symmetry group
    Orbit(Element generating_element,
          const SymGroup& generating_group, 
          const SymCompare<Element>& sym_compare);
    
    /// \brief Construct an Orbit from a generating_element Element, using provided symmetry op
    template<typename SympOpIterator>
    Orbit(Element generating_element, 
          SympOpIterator begin,
          SympOpIterator end, 
          const SymCompare<Element>& sym_compare);
    
    const_iterator begin() const {
      return m_element.cbegin();
    }
    
    const_iterator end() const {
      return m_element.cend();
    }
    
    const_iterator cbegin() const {
      return m_element.cbegin();
    }
    
    const_iterator cend() const {
      return m_element.cend();
    }
    
    size_type size() const {
      return m_element.size();
    }
    
    /// \brief Identical to element(0)
    const Element& prototype() const {
      return m_element[0];
    }
    
    /// \brief Return Element at index, without bounds checking
    const Element& operator[](size_type index) const {
      return element(index);
    }
    
    /// \brief Equivalent to operator[](size_type index) const
    const Element& element(size_type index) const {
      return m_element[index];
    }
    
    /// \brief const Access vector of Element
    const std::vector<Element>& elements() const {
      return m_element;
    }
    
    /// \brief Return the equivalence map for element[index]
    ///
    /// \returns a pair of const_iterators, begin and end, over SymOp such that
    /// element(index) compares equivalent to op*prototype()
    std::pair<const_iterator, const_iterator> equivalence_map(size_type index) const {
      return std::make_pair(m_equivalence_map.begin(), m_equivalence_map.end());
    }
    
    /// \brief Return the SymCompare functor reference
    ///
    /// - implements symmetry properties of this orbit
    const SymCompare& sym_compare() const {
      return *m_sym_compare;
    }
    
    /// \brief Apply symmetry to Orbit
    Orbit& apply_sym(const SymOp &op);
    
    private:
    
    /// \brief Construct an Orbit from a generating_element Element, using provided symmetry rep
    template<typename SymOpIterator>
    void _construct(Element generating_element, 
                    SymOpIterator begin,
                    SymOpIterator end);
    
    /// \brief All symmetrically equivalent elements (excluding translations)
    std::vector<Element> m_element;
    
    /// \brief element(i) is equivalent to m_equivalence_map(i, j)*prototype() for all j
    std::vector<std::vector<SymOp> > m_equivalence_map;
    
    /// \brief Functor used to check compare Element, including symmetry rules,
    /// and make canonical forms
    notstd::cloneable_ptr<SymCompare<Element> > m_sym_compare;
    
  };
  
  
  /* -- Orbit Definitions ------------------------------------- */
  
  /// \brief Construct an Orbit from a generating_element Element, using provided Group
  /// 
  /// \param generating_element The element used to generate equivalents
  /// \param generating_group The group used for generating the orbit
  /// \param sym_compare Binary functor that implements symmetry properties
  ///
  template<typename Element>
  Orbit<Element>::Orbit(Element generating_element, 
                        const SymGroup& generating_group, 
                        const SymCompare<Element>& sym_compare) :
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
                        const SymCompare<Element>& sym_compare) :
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
  Orbit<Element>::_construct(Element generating_element, 
                             SympOpIterator begin,
                             SympOpIterator end) {
    
    // add all op*generating_element
    for(auto it=begin; it!=end; ++it) {
      m_element.push_back(m_sym_compare->prepare(copy_apply(*it, generating_element)));
    }
    
    // sort
    std::sort(m_element.begin(), m_element.end(), *m_sym_compare);
    
    // keep uniques
    auto last = std::unique(m_element.begin(), m_element.end());
    m_element.erase(last, m_element.end());
    
    // create equivalence map
    for(auto it=begin; it!=end; ++it) {
      auto equiv = m_sym_compare->prepare(copy_apply(*it, prototype()));
      equivalence_map[find_index(m_element, equiv)].push_back(*it);
    }
  }
  
  /// \brief Apply symmetry to Orbit
  template<typename Element>
  Orbit<Element>& Orbit<Element>::apply_sym(const SymOp &op) {
    
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