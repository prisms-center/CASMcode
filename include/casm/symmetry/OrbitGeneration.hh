#ifndef CASM_OrbitGeneration
#define CASM_OrbitGeneration

#include <utility>
#include <set>

#include "casm/symmetry/SymGroup.hh"

namespace CASM {

  class PrimClex;
  class SymGroup;

  /** \defgroup OrbitGeneration

      \brief Helpers for generating Orbit
  */

  template<typename _OrbitType>
  struct OrbitGeneratorCompare;

  template<typename _OrbitType>
  struct CanonicalGenerator;

  /// \brief An std::set of Orbit
  ///
  /// \ingroup OrbitGeneration
  ///
  template<typename OrbitType>
  using OrbitGeneratorSet = std::set<typename OrbitType::Element, OrbitGeneratorCompare<OrbitType> >;

  /// \brief Data structure that holds canonical generating elements and
  /// can then make sorted orbits
  ///
  /// \ingroup OrbitGeneration
  ///
  template<typename _OrbitType>
  struct OrbitGenerators {

    typedef _OrbitType OrbitType;
    typedef typename OrbitType::Element Element;
    typedef typename OrbitType::SymCompareType SymCompareType;

    OrbitGenerators(const SymGroup &_group,
                    const SymCompareType &_sym_compare);

    /// \brief Try inserting an element, after generating the canonical form
    std::pair<typename OrbitGeneratorSet<OrbitType>::iterator, bool> insert(const Element &test);

    /// \brief Try inserting an element, assuming it is in canonical form
    std::pair<typename OrbitGeneratorSet<OrbitType>::iterator, bool> insert_canonical(const Element &test);

    /// \brief Construct Orbit from all generating elements
    template<typename OrbitOutputIterator>
    OrbitOutputIterator make_orbits(OrbitOutputIterator result);

    /// \brief Construct Orbit from all generating elements, including PrimClex pointer
    template<typename OrbitOutputIterator>
    OrbitOutputIterator make_orbits(OrbitOutputIterator result, const PrimClex &primclex);

    const SymGroup &group;
    const SymCompareType &sym_compare;

  private:

    OrbitGeneratorCompare<OrbitType> m_element_compare;
    CanonicalGenerator<OrbitType> m_generate_canonical;

  public:

    OrbitGeneratorSet<OrbitType> elements;

  };

  /// \brief Compare concept functor for canonical generating elements
  ///
  /// - Uses SymCompareType::compare
  ///
  /// \ingroup OrbitGeneration
  ///
  template<typename _OrbitType>
  struct OrbitGeneratorCompare {

    typedef _OrbitType OrbitType;
    typedef typename OrbitType::Element Element;
    typedef typename OrbitType::SymCompareType SymCompareType;

    const SymCompareType &sym_compare;

    OrbitGeneratorCompare(const SymCompareType &_sym_compare);

    bool operator()(const Element &A, const Element &B) const;
  };

  /// \brief Functor to find the canonical generating element for an orbit
  ///
  /// - Uses generating SymGroup, SymCompareType::prepare, SymCompareType::compare
  ///
  /// \ingroup OrbitGeneration
  ///
  template<typename _OrbitType>
  struct CanonicalGenerator {

    typedef _OrbitType OrbitType;
    typedef typename OrbitType::Element Element;
    typedef typename OrbitType::SymCompareType SymCompareType;

    const SymGroup &generating_group;
    const SymCompareType &sym_compare;

    CanonicalGenerator(
      const SymGroup &_generating_group,
      const SymCompareType &_sym_compare);

    /// \brief Applies symmetry to return an equivalent Element in a canonical form
    Element operator()(const Element &e) const;

    /// \brief After using call operator, this can be checked
    const SymOp &to_canonical() const;

    /// \brief After using call operator, this can be checked
    const SymOp &from_canonical() const;

  private:

    mutable const SymOp *m_to_canonical;
  };

  /// \brief Functor to find to check if element is in canonical form
  ///
  /// - Uses generating SymGroup, SymCompareType::prepare, SymCompareType::compare
  ///
  /// \ingroup OrbitGeneration
  ///
  template<typename _OrbitType>
  struct IsCanonical {

    typedef _OrbitType OrbitType;
    typedef typename OrbitType::Element Element;
    typedef typename OrbitType::SymCompareType SymCompareType;

    const SymGroup &generating_group;
    const SymCompareType &sym_compare;

    IsCanonical(
      const SymGroup &_generating_group,
      const SymCompareType &_sym_compare);

    /// \brief Applies symmetry to check if any Element is greater than e
    bool operator()(const Element &e) const;

  };

}

#endif
