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


  /// \brief CRTP base class for implementing element and orbit comparison
  ///
  /// Derived needs to implement the following *private* methods:
  /// - Element Derived::spatial_prepare_impl(const Element &A) const
  ///   - position/orient Element in a canonical way for comparison
  ///     (i.e., translate cluster so first site is at [0,0,0] cell)
  /// - Element Derived::representation_prepare_impl(const Element &A) const
  ///   - construct one canonical representation of Element amongst all equivalent representations
  ///     (i.e. sort sites of a cluster), used for comparison
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
  /// - void Derived::copy_apply_sym(SymOp const &op, Element obj) const
  ///   - defines application of SymOp to object of type Element
  ///   - default uses standalone copy_apply() function
  ///   - allows flexibility for working outside the typical MasterSymGroup/SymGroupRep framework
  ///     or when symmetry application is ambiguous (like applying symmetry to a Eigen::VectorXd)
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

    using MostDerived = typename Base::MostDerived;
    using Base::derived;
    using Element = typename traits<MostDerived>::Element;

    /// \brief Applies SymOp to Element. Default performs standalone copy_apply
    Element copy_apply(SymOp const &op, Element obj) const {
      return derived().copy_apply_impl(op, obj);
    }

    /// \brief Prepare an element for comparison via an isometric affine transformation
    ///
    /// - For instance, translate a cluster so comparison may be performed more efficiently.
    /// - The transformation that is applied is stored and available as `spatial_transform`
    /// - The `spatial_transform` is also included in the symmetry operations stored in the
    ///   Orbit equivalence_map (i.e. equivalence_map_element = spatial_transform * generating_group_op)
    /// - Returns pair such that pair.first = apply_sym(pair.second, obj)
    ///
    Element spatial_prepare(Element obj) const {
      return derived().spatial_prepare_impl(obj);
    }

    /// \brief Prepare an element for comparison via transformation of its internal representation
    ///
    /// - For instance, sort sites of a cluster so comparison may be
    ///   performed more efficiently
    ///
    Element representation_prepare(Element obj) const {
      return derived().representation_prepare_impl(obj);
    }

    /// \brief Prepare an element for comparison via representation_prepare(), followed by spatial_prepare()
    ///
    Element prepare(Element obj) const {
      return spatial_prepare(representation_prepare(obj));
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
    bool invariants_compare(const Element &A, const Element &B) const {
      return derived().invariants_compare_impl(A, B);
    }

    /// \brief Check equivalence of invariants
    ///
    /// \returns \code !invariants_compare(A,B) && !invariants_compare(B,A) \endcode
    bool invariants_equal(const Element &A, const Element &B) const {
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

    /// \brief Access spatial transform that was used during most recent spatial preparation of an element
    SymOp const &spatial_transform() const {
      return m_spatial_transform;
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
      if(this->invariants_compare(A, B)) {
        return true;
      }
      if(this->invariants_compare(B, A)) {
        return false;
      }

      // next compare A and B
      return this->compare(A, B);
    }

    /// \brief Default spatial_prepare does nothing
    Element spatial_prepare_impl(Element obj) const {
      return obj;
    }

    /// \brief Applies SymOp to Element. Default performs standalone copy_apply
    Element copy_apply_impl(SymOp const &op, Element obj) const {
      return CASM::copy_apply(op, obj);
    }

    /// \brief Transforms the SymCompare object, default does nothing
    SymCompare &apply_sym_impl(const SymOp &op) {
      return *this;
    }

    /// \brief Spatial transform that reproduces most recent application of SymCompare::spatial_prepare()
    /// - Default SymOp constructor initializes to identity
    mutable SymOp m_spatial_transform;
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
