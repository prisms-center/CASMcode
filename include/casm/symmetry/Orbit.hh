#ifndef CASM_Orbit
#define CASM_Orbit

#include <vector>

#include <boost/iterator/transform_iterator.hpp>

#include "casm/misc/CASM_math.hh"
#include "casm/misc/Comparisons.hh"
#include "casm/container/multivector.hh"
#include "casm/symmetry/OrbitDecl.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/SymGroup.hh"

namespace CASM {
  class SymGroup;

  // -- Orbit -------------------------------------

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
  /// The following relationships will be valid:
  ///    element(i) = sym_compare.copy_apply(equivalence_map()[i][j], prototype()), for all j < equivalence_map()[i].size()
  ///    equivalence_map()[i][j] = t * g;
  ///      where g is a generating group element,
  ///      and t is the "spatial_transform" defined by the SymCompareType for equivalent elements such that:
  ///    sym_compare.representation_prepare(element(i)) ==
  ///      copy_apply(t, sym_compare.representation_prepare(copy_apply(g, prototype)))
  ///
  /// \ingroup Clusterography
  ///
  template<typename _SymCompareType>
  class Orbit : public Comparisons<CRTPBase<Orbit<_SymCompareType>>> {

  public:

    using size_type = unsigned int;
    using Element = typename _SymCompareType::Element;
    using InvariantsType = typename _SymCompareType::InvariantsType;
    using SymCompareType = _SymCompareType;
    using const_iterator = typename std::vector<Element>::const_iterator;
    using const_symop_iterator = typename std::vector<SymOp>::const_iterator;

    /// \brief Construct an Orbit from a generating_element Element, using provided symmetry group
    Orbit(Element generating_element,
          SymGroup const &generating_group,
          SymCompareType const &sym_compare);


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
    ///
    /// - May not be prepared
    const Element &operator[](size_type index) const {
      return element(index);
    }

    /// \brief Equivalent to operator[](size_type index) const
    ///
    /// - May not be prepared
    const Element &element(size_type index) const {
      return m_element[index];
    }

    /// \brief const Access vector of Element
    ///
    /// - May not be prepared
    const std::vector<Element> &elements() const {
      return m_element;
    }

    /// \brief Return the equivalence map
    ///
    /// \returns element(i) compares equivalent to prototype().copy_apply(equivalence_map[i][j]) for all j
    ///
    const multivector<SymOp>::X<2> &equivalence_map() const {
      return m_equivalence_map;
    }

    /// \brief Return the equivalence map for element[index]
    ///
    /// \returns a pair of const_iterators, begin and end, over SymOp such that
    /// element(index) compares equivalent to prototype().copy_apply(op)
    std::pair<const_symop_iterator, const_symop_iterator> equivalence_map(size_type index) const {
      return std::make_pair(m_equivalence_map[index].begin(), m_equivalence_map[index].end());
    }

    /// \brief Return the canonization symmetry representation ID
    SymGroupRepID canonization_rep_ID() const {
      if(m_canonization_rep_ID.empty())
        _construct_canonization_rep();
      return m_canonization_rep_ID;
    }

    /// \brief Find element in Orbit
    ///
    /// - Assumes 'e' is 'prepared', uses SymCompare<Element>::intra_orbit_equal
    ///   to check equivalence
    const_iterator find(const Element &e) const {
      return std::find_if(begin(), end(), [&](const Element & B) {
        return m_sym_compare.equal(e, m_sym_compare.prepare(B));
      });
    }

    /// \brief Check if element is in Orbit
    ///
    /// - Assumes 'e' is 'prepared', uses SymCompare<Element>::intra_orbit_equal
    ///   to check equivalence
    bool contains(const Element &e) const {
      return this->find(e) != end();
    }

    /// \brief Return the generating SymGroup
    const SymGroup &generating_group() const {
      return m_generating_group;
    }

    /// \brief Return the SymCompare functor reference
    ///
    /// - implements symmetry properties of this orbit
    const SymCompareType &sym_compare() const {
      return m_sym_compare;
    }

    InvariantsType const &invariants() const {
      return m_invariants;
    }

    /// \brief Apply symmetry to Orbit
    Orbit &apply_sym(const SymOp &op);

    /// \brief Compare orbits, using SymCompareType::inter_orbit_compare
    bool operator<(const Orbit &B) const;

  private:

    void _construct_canonization_rep() const;

    /// \brief Construct an Orbit from a generating_element Element, using provided symmetry rep
    template<typename SymOpIterator>
    void _construct(Element generating_element,
                    SymOpIterator begin,
                    SymOpIterator end);

    /// \brief All symmetrically equivalent elements (excluding those that SymCompare equivalent)
    ///
    std::vector<Element> m_element;

    /// \brief element(i) compares equivalent to prototype().copy_apply(m_equivalence_map[i][j]) for all j
    multivector<SymOp>::X<2> m_equivalence_map;

    /// \brief Group used to generate the orbit
    SymGroup m_generating_group;

    /// \brief ID of symmetry representation that describes the effect of each SymOp with respect to the canonical equivalent
    mutable SymGroupRepID m_canonization_rep_ID;

    /// \brief Functor used to check compare Element, including symmetry rules,
    /// and make canonical forms
    SymCompareType m_sym_compare;

    /// \brief Orbit invariants
    InvariantsType m_invariants;
  };


  /// \brief Iterator over Generators (potential prototypes) and insert resulting orbits into 'result' iterator
  template<typename GeneratorIterator, typename SymCompareType, typename OrbitOutputIterator>
  OrbitOutputIterator make_orbits(
    GeneratorIterator gen_begin,
    GeneratorIterator gen_end,
    const SymGroup &generating_group,
    const SymCompareType &sym_compare,
    OrbitOutputIterator result) {

    using OrbitType = typename OrbitOutputIterator::container_type::value_type;
    for(; gen_begin != gen_end; ++gen_begin) {
      *(result++) = OrbitType(*gen_begin, generating_group, sym_compare);
    }
    return result;
  }

  // -- Orbit Helpers --------------------

  /// \brief Find orbit containing an element in a range of Orbit
  template<typename OrbitIterator, typename Element>
  OrbitIterator find_orbit(OrbitIterator begin, OrbitIterator end, Element e);

  struct GetPrototype {

    template<typename OrbitType>
    typename OrbitType::Element const &operator()(const OrbitType &orbit) const {
      return orbit.prototype();
    }
  };

  template<typename OrbitIterator>
  using PrototypeIterator = boost::transform_iterator<GetPrototype, OrbitIterator>;

  /// Convert an Orbit iterator to a prototype iterator
  template<typename OrbitIterator>
  PrototypeIterator<OrbitIterator> prototype_iterator(OrbitIterator orbit_it) {
    return boost::make_transform_iterator(orbit_it, GetPrototype());
  }


  struct GetInvariants {

    template<typename OrbitType>
    typename OrbitType::InvariantsType const &operator()(const OrbitType &orbit) const {
      return orbit.invariants();
    }
  };

  template<typename OrbitIterator>
  using InvariantsIterator = boost::transform_iterator<GetInvariants, OrbitIterator>;

  /// Convert an Orbit iterator to an invariants iterator
  template<typename OrbitIterator>
  InvariantsIterator<OrbitIterator> invariants_iterator(OrbitIterator orbit_it) {
    return boost::make_transform_iterator(orbit_it, GetInvariants());
  }
}

#endif
