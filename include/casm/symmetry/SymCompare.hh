#ifndef CASM_SymCompare
#define CASM_SymCompare

#include <utility>
#include <memory>
#include "casm/CASM_global_definitions.hh"
#include "casm/symmetry/SymGroup.hh"

namespace CASM {

  /* -- SymCompare Declarations --------------------------- */

  class SymOp;

  /// \brief CRTP base class for implementing element and orbit comparison
  ///
  /// Derived needs to implement the following *private* methods:
  /// - Element Derived::prepare_impl(const Element &A) const
  ///   - use to put Element into a canonical comparison form, i.e. sort a cluster
  /// - bool Derived::compare_impl(const Element &A, const Element &B) const
  ///   - used to identify unique equivalents
  /// - bool Derived::invariants_compare_impl(const Element &A, const Element &B) const
  ///   - first order comparison of orbits
  ///   - speeds up identifying orbits that might contain an element
  /// - bool Derived::inter_orbit_compare_impl(const Element &A, const Element &B) const
  ///   - breaks 'invariants_compare' ties
  /// - void Derived::apply_sym_impl(const SymOp &op) const
  ///   - not sure if this is needed... use to apply_sym to the SymCompare object itself
  ///
  /// The ClusterSymCompare hierarchy:
  /// - SymCompare
  ///   - ClusterSymCompare (implements 'invariants_compare_impl', 'inter_orbit_compare_impl', and 'apply_sym_impl')
  ///     - IntegralClusterSymCompare (implements 'compare_impl')
  ///       - LocalSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///       - PrimPeriodicSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///       - ScelPeriodicSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///
  template<typename Derived>
  class SymCompare {

  public:

    typedef typename CASM_TMP::traits<Derived>::MostDerived MostDerived;
    typedef typename CASM_TMP::traits<Derived>::Element Element;
    typedef typename CASM_TMP::traits<Derived>::InvariantsType InvariantsType;

    /// \brief Prepare an element for comparison
    ///
    /// - For instance, sort and/or translate a cluster so comparison may be
    ///   performed more efficiently
    Element prepare(Element obj) const {
      return derived().prepare_impl(obj);
    }

    /// \brief Orders 'prepared' elements
    ///
    /// - Returns 'true' to indicate A < B
    /// - Assumes elements are 'prepared' before being compared
    bool compare(const Element &A, const Element &B) const {
      return derived().compare_impl(A, B);
    }

    /// \brief Check equivalence of 'prepared' elements
    ///
    /// \returns \code !compare(A,B) && !compare(B,A) \endcode
    ///
    /// - Assumes elements are 'prepared' before being compared
    bool equal(const Element &A, const Element &B) const {
      return !compare(A, B) && !compare(B, A);
    }

    /// \brief Defines an order for elements that have different invariants
    ///
    /// - Returns 'true' to indicate A < B
    /// - Equivalence is indicated by \code !invariants_compare(A,B) && !invariants_compare(B,A) \endcode
    /// - Assumes elements are 'prepared' before being compared
    bool invariants_compare(const InvariantsType &A, const InvariantsType &B) const {
      return derived().invariants_compare_impl(A, B);
    }

    /// \brief Check equivalence of invariants
    ///
    /// \returns \code !invariants_compare(A,B) && !invariants_compare(B,A) \endcode
    bool invariants_equal(const InvariantsType &A, const InvariantsType &B) const {
      return !invariants_compare(A, B) && !invariants_compare(B, A);
    }

    /// \brief Orders orbit prototypes, breaking invariants_compare ties
    ///
    /// - Returns 'true' to indicate A < B
    /// - Equivalence is indicated by \code !compare(A,B) && !compare(B,A) \endcode
    /// - Assumes elements are in canonical form
    bool inter_orbit_compare(const Element &A, const Element &B) const {
      return derived().inter_orbit_compare_impl(A, B);
    }

    /// \brief Check equivalence of prototypes in different orbit
    ///
    /// \returns \code !compare(A,B) && !compare(B,A) \endcode
    ///
    /// - Assumes elements are in canonical form
    bool inter_orbit_equal(const Element &A, const Element &B) const {
      return !inter_orbit_compare(A, B) && !inter_orbit_compare(B, A);
    }

    /// \brief Apply symmetry to this
    ///
    /// - For instance, this might change the direction of periodicity, affecting
    ///   the outcome of 'prepare' and/or 'compare'
    SymCompare &apply_sym(const SymOp &op) {
      derived().apply_sym_impl(op);
      return *this;
    }

  protected:

    SymCompare() {}

    MostDerived &derived() {
      return *static_cast<MostDerived *>(this);
    }

    const MostDerived &derived() const {
      return *static_cast<const MostDerived *>(this);
    }
  };

  /// \brief Template class to be specialized for comparisons with aperiodic symmetry
  template<typename Element>
  class LocalSymCompare {};

  /// \brief Template class to be specialized for comparisons with periodic symmetry
  /// of the primitive lattice
  template<typename Element>
  class PrimPeriodicSymCompare {};

  /// \brief Template class to be specialized for comparisons with periodic symmetry
  /// of the supercell lattice
  template<typename Element>
  class ScelPeriodicSymCompare {};

  /// \brief Return subgroup that leaves an element unchanged
  ///
  /// All SymOp such that:
  /// \code
  /// Element e = sym_compare.prepare(generating_element);
  /// Element test = sym_compare.prepare(copy_apply(op, e));
  /// sym_compare.equal(e, test) == true
  /// \endcode
  template<typename Element, typename SymCompareType>
  SymGroup invariant_subgroup(const Element &element,
                              const SymGroup &generating_grp,
                              const SymCompareType &sym_compare) {
    Element e(sym_compare.prepare(element));
    SymGroup result;
    for(const auto &op : generating_grp) {
      if(sym_compare.equal(e, sym_compare.prepare(copy_apply(op, e)))) {
        result.push_back(op);
      }
    }
    return result;
  }

}

#endif
