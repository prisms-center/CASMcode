#ifndef CASM_Orbit
#define CASM_Orbit

#include <vector>

#include "casm/misc/cloneable_ptr.hh"
#include "casm/symmetry/SymCompare.hh"

namespace CASM {

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
          const SymGroup &generating_group,
          const SymCompare<Element> &sym_compare);

    /// \brief Construct an Orbit from a generating_element Element, using provided symmetry op
    template<typename SympOpIterator>
    Orbit(Element generating_element,
          SympOpIterator begin,
          SympOpIterator end,
          const SymCompare<Element> &sym_compare);

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
    const Element &prototype() const {
      return m_element[0];
    }

    /// \brief Return Element at index, without bounds checking
    const Element &operator[](size_type index) const {
      return element(index);
    }

    /// \brief Equivalent to operator[](size_type index) const
    const Element &element(size_type index) const {
      return m_element[index];
    }

    /// \brief const Access vector of Element
    const std::vector<Element> &elements() const {
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
    const SymCompare &sym_compare() const {
      return *m_sym_compare;
    }

    /// \brief Apply symmetry to Orbit
    Orbit &apply_sym(const SymOp &op);

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
}

#endif
