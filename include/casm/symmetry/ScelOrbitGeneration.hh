#ifndef CASM_ScelOrbitGeneration
#define CASM_ScelOrbitGeneration

// TODO: this does not belong in casm/symmetry

#include "casm/clusterography/ClusterSymCompareDecl.hh"
#include "casm/symmetry/OrbitDecl.hh"
#include "casm/symmetry/PermuteIterator.hh"

namespace CASM {

  class Supercell;

  /** \defgroup OrbitGeneration

      \brief Helpers for generating Orbit
  */


  /// \brief Functor to find the canonical generating element for an orbit in a Supercell
  ///
  /// - Uses Supercell permute group and crystallography_tol
  ///
  /// \ingroup OrbitGeneration
  ///
  template<typename _ElementType>
  class ScelCanonicalGenerator {
  public:

    typedef _ElementType Element;
    typedef ScelPeriodicSymCompare<Element> SymCompareType;
    typedef Orbit<SymCompareType> OrbitType;

    ScelCanonicalGenerator(const Supercell &_scel);

    const Supercell &supercell() const;

    const SymCompareType &sym_compare() const;

    /// \brief Applies symmetry to return an equivalent Element in a canonical form
    ///
    /// - Use [supercell().permute_begin(), supercell().permute_end()) to canonicalize
    Element operator()(const Element &e) const;

    /// \brief Applies symmetry to return an equivalent Element in a canonical form
    ///
    /// - For use with a container of PermuteIterator
    /// - Use [begin, end) to canonicalize
    template<typename PermuteIteratorIt>
    Element operator()(const Element &e, PermuteIteratorIt begin, PermuteIteratorIt end) const;

    /// \brief After using call operator, this can be checked
    PermuteIterator to_canonical() const;

    /// \brief After using call operator, this can be checked
    PermuteIterator from_canonical() const;

  private:

    const Supercell *m_scel;
    SymCompareType m_sym_compare;
    mutable PermuteIterator m_to_canonical;
  };

  /// \brief Functor to find to check if element is in canonical form
  ///
  /// - Uses generating SymGroup, SymCompareType::prepare, SymCompareType::compare
  ///
  /// \ingroup OrbitGeneration
  ///
  template<typename _ElementType>
  struct ScelIsCanonical {

    typedef _ElementType Element;
    typedef ScelPeriodicSymCompare<Element> SymCompareType;

    ScelIsCanonical(const Supercell &_scel);

    const Supercell &supercell() const;

    const SymCompareType &sym_compare() const;

    /// \brief Applies symmetry to check if any Element is greater than e
    /// - Use [supercell().permute_begin(), supercell().permute_end()) to canonicalize
    bool operator()(const Element &e) const;

    /// \brief Applies symmetry to check if any Element is greater than e
    template<typename PermuteIteratorIt>
    bool operator()(
      const Element &e,
      PermuteIteratorIt begin,
      PermuteIteratorIt end) const;

    const Supercell *m_scel;
    SymCompareType m_sym_compare;
  };

}

#endif
