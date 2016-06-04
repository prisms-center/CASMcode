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
  /// map giving symmetry operations that map the prototype onto the
  /// equivalents.
  ///
  /// Element and orbit comparison is done via a SymCompareType functor, which
  /// includes any necessary tolerance for floating point comparison. See `SymCompare`
  /// for how to implement the necessary methods.
  ///
  /// \ingroup Clusterography
  ///
  template<typename _Element, typename _SymCompareType = PrimPeriodicSymCompare<_Element> >
  class Orbit : public Comparisons<Orbit<_Element, _SymCompareType> > {

  public:

    typedef unsigned int size_type;
    typedef _Element Element;
    typedef _SymCompareType SymCompareType;
    typedef typename std::vector<Element>::const_iterator const_iterator;

    Orbit() {}

    /// \brief Construct an Orbit from a generating_element Element, using provided symmetry group
    Orbit(Element generating_element,
          const SymGroup &generating_group,
          const SymCompareType &sym_compare);

    /// \brief Construct an Orbit from a generating_element Element, using provided symmetry op
    template<typename SympOpIterator>
    Orbit(Element generating_element,
          SympOpIterator begin,
          SympOpIterator end,
          const SymCompareType &sym_compare);

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

    /// \brief Find element in Orbit
    ///
    /// - Assumes 'e' is 'prepared', uses SymCompare<Element>::intra_orbit_equal
    ///   to check equivalence
    const_iterator find(const Element &e) {
      return std::find_if(begin(), end(), [&](const Element & B) {
        return m_sym_compare->intra_orbit_equal(e, B);
      });
    }

    /// \brief Check if element is in Orbit
    ///
    /// - Assumes 'e' is 'prepared', uses SymCompare<Element>::intra_orbit_equal
    ///   to check equivalence
    bool contains(const Element &e) {
      return this->find(e) != end();
    }

    /// \brief Return the SymCompare functor reference
    ///
    /// - implements symmetry properties of this orbit
    const SymCompareType &sym_compare() const {
      return m_sym_compare;
    }

    /// \brief Apply symmetry to Orbit
    Orbit &apply_sym(const SymOp &op);

    /// \brief Compare orbits, using SymCompareType::inter_orbit_compare
    bool operator<(const Orbit &B) const {
      return m_sym_compare.inter_orbit_compare(prototype(), B.prototype());
    }

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
    SymCompareType m_sym_compare;

  };


  /// \brief Find orbit containing an element in a range of Orbit
  template<typename OrbitIterator, typename Element>
  OrbitIterator find_orbit(OrbitIterator begin, OrbitIterator end, Element e);
}

#endif
