#ifndef CASM_SymCompare
#define CASM_SymCompare

#include <utility>
#include <memory>
#include "casm/CASM_global_definitions.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/OrbitDecl.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

  /* -- SymCompare Declarations --------------------------- */

  class SymOp;
  class UnitCell;
  class Coordinate;

  /// \brief Traits class for AperiodicSymCompare
  template<typename _Element>
  struct traits<AperiodicSymCompare<_Element>> {
    typedef _Element Element;
    typedef AperiodicSymCompare<Element> MostDerived;
  };

  /// \brief Traits class for PrimPeriodicSymCompare
  template<typename _Element>
  struct traits<PrimPeriodicSymCompare<_Element>> {
    typedef _Element Element;
    typedef PrimPeriodicSymCompare<Element> MostDerived;
  };

  /// \brief Traits class for ScelPeriodicSymCompare
  template<typename _Element>
  struct traits<ScelPeriodicSymCompare<_Element>> {
    typedef _Element Element;
    typedef ScelPeriodicSymCompare<Element> MostDerived;
  };

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
  ///
  /// Derived should optionally implement:
  /// - bool Derived::inter_orbit_compare_impl(const Element &A, const Element &B) const
  ///   - breaks 'invariants_compare' ties
  ///   - default uses this->compare(A, B)
  /// - void Derived::apply_sym_impl(const SymOp &op) const
  ///   - not sure if this is needed... use to apply_sym to the SymCompare object itself
  ///   - default does nothing
  ///
  /// The SymCompare hierarchy for ClusterTypes:
  /// - SymCompare
  ///   - ClusterSymCompare (implements 'invariants_compare_impl', 'inter_orbit_compare_impl', and 'apply_sym_impl')
  ///     - LocalSymCompare<ClusterType> (implements 'compare_impl')
  ///       - LocalSymCompare<ClusterType> (implements 'prepare_impl')
  ///       - PrimPeriodicSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///       - ScelPeriodicSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///
  template<typename Base>
  class SymCompare : public Base {

  public:

    typedef typename Base::MostDerived MostDerived;
    using Base::derived;
    typedef typename traits<MostDerived>::Element Element;
    typedef typename traits<Element>::InvariantsType InvariantsType;


    SymCompare() : m_integral_tau(Eigen::Vector3l::Zero(3)) {}

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
    /// - Equivalence is indicated by \code !inter_orbit_compare(A,B) && !inter_orbit_compare(B,A) \endcode
    /// - Assumes elements are in canonical form
    bool inter_orbit_compare(const Element &A, const Element &B) const {
      return derived().inter_orbit_compare_impl(A, B);
    }

    /// \brief Check equivalence of prototypes in different orbit
    ///
    /// \returns \code !inter_orbit_compare(A,B) && !inter_orbit_compare(B,A) \endcode
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

    /// \brief Access transform that took cluster from unprepared to prepared state
    std::unique_ptr<SymOpRepresentation> canonical_transform(Element const  &obj) const {
      return derived().canonical_transform_impl(obj);
    }

    /// \brief Access integral adjustment shift due to varying symmetry of object vs. generating group
    const UnitCell integral_tau() const {
      return m_integral_tau;
    }

    /// \brief Access SymOp adjustment due to varying symmetry of object vs. generating group
    const SymOp translation(const Structure &prim) const {
      Coordinate tau(m_integral_tau.cast<double>(), prim.lattice(), FRAC);
      return SymOp::translation(tau.const_cart());
    }

  protected:

    /// \brief Orders orbit prototypes, breaking invariants_compare ties
    ///
    /// - Returns 'true' to indicate A < B
    /// - Equivalence is indicated by \code !compare(A,B) && !compare(B,A) \endcode
    ///
    /// Implementation:
    /// - First, check invariants_compare
    /// - Break ties with, compare(A, B)
    bool inter_orbit_compare_impl(const Element &A, const Element &B) const {

      // first compare invariants
      if(this->invariants_compare(A.invariants(), B.invariants())) {
        return true;
      }
      if(this->invariants_compare(B.invariants(), A.invariants())) {
        return false;
      }

      // next compare A and B
      return this->compare(A, B);
    }

    /// \brief Transforms the SymCompare object, default does nothing
    SymCompare &apply_sym_impl(const SymOp &op) {
      return *this;
    }

    mutable UnitCell m_integral_tau;
  };


  /// \brief CRTP Base class for types that should be SymComparable
  template<typename _Base>
  class SymComparable : public _Base {
  public:

    typedef _Base Base;
    typedef typename Base::MostDerived MostDerived;
    using Base::derived;
    typedef typename traits<MostDerived>::InvariantsType InvariantsType;

    const InvariantsType &invariants() const {
      if(!m_invariants) {
        m_invariants = notstd::make_cloneable<InvariantsType>(derived());
      }
      return *m_invariants;
    }

  protected:

    void reset_invariants() {
      m_invariants.reset();
    }

    mutable notstd::cloneable_ptr<InvariantsType> m_invariants;
  };

}

#endif
