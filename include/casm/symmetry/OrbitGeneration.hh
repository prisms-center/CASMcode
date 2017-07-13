#ifndef CASM_OrbitGeneration
#define CASM_OrbitGeneration

#include <utility>
#include <set>

namespace CASM {

  template<typename _Element, typename _SymCompareType>
  class Orbit;

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
                    const SymCompareType &_sym_compare) :
      group(_group),
      sym_compare(_sym_compare),
      m_element_compare(sym_compare),
      m_generate_canonical(group, sym_compare),
      elements(m_element_compare) {
    }

    /// \brief Try inserting an element, after generating the canonical form
    std::pair<typename OrbitGeneratorSet<OrbitType>::iterator, bool> insert(const Element &test) {
      return elements.insert(m_generate_canonical(test));
    }

    /// \brief Try inserting an element, assuming it is in canonical form
    std::pair<typename OrbitGeneratorSet<OrbitType>::iterator, bool> insert_canonical(const Element &test) {
      return elements.insert(test);
    }

    /// \brief Construct Orbit from all generating elements
    template<typename OrbitOutputIterator>
    OrbitOutputIterator make_orbits(OrbitOutputIterator result) {

      // orbit generating elements should already by sorted
      std::vector<OrbitType> orbits;
      for(const auto &e : elements) {
        orbits.emplace_back(e, group, sym_compare);
      }
      // output Orbits
      return std::move(orbits.begin(), orbits.end(), result);
    }

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

    OrbitGeneratorCompare(const SymCompareType &_sym_compare) :
      sym_compare(_sym_compare) {}

    bool operator()(const Element &A, const Element &B) const {
      return sym_compare.inter_orbit_compare(A, B);
    };
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
      const SymCompareType &_sym_compare) :
      generating_group(_generating_group),
      sym_compare(_sym_compare),
      m_to_canonical(nullptr) {}

    /// \brief Applies symmetry to return an equivalent Element in a canonical form
    Element operator()(const Element &e) const {
      Element result = sym_compare.prepare(e);
      for(const auto &op : generating_group) {
        auto test = sym_compare.prepare(copy_apply(op, e));
        if(sym_compare.compare(result, test)) {
          result = test;
          m_to_canonical = &op;
        }
      }
      return result;
    }

    /// \brief After using call operator, this can be checked
    const SymOp &to_canonical() const {
      return *m_to_canonical;
    }

    /// \brief After using call operator, this can be checked
    const SymOp &from_canonical() const {
      return to_canonical().inverse();
    }

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
      const SymCompareType &_sym_compare) :
      generating_group(_generating_group),
      sym_compare(_sym_compare) {}

    /// \brief Applies symmetry to check if any Element is greater than e
    bool operator()(const Element &e) const {
      return std::none_of(
               generating_group.begin(),
               generating_group.end(),
      [&](const SymOp & op) {
        auto test = sym_compare.prepare(copy_apply(op, e));
        return sym_compare.compare(e, test);
      });
    }

  };

}

#endif
