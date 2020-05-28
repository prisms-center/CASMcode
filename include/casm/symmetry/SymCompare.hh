#ifndef CASM_SymCompare
#define CASM_SymCompare

#include <utility>
#include <memory>
#include "casm/global/definitions.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

  class SymOp;
  class SymOpRepresentation;

  /// \brief CRTP base class for implementing element and orbit comparison
  ///
  /// Derived needs to implement the following *private* methods:
  /// - Element Derived::spatial_prepare_impl(const Element &A) const
  ///   - position/orient Element in a canonical way for comparison
  ///     (i.e., translate cluster so first site is at [0,0,0] cell)
  /// - SymOp const &Derived::spatial_transform_impl() const;
  ///   - Access spatial transform that was used during most recent spatial preparation of an
  ///     element
  /// - Element Derived::representation_prepare_impl(const Element &A) const
  ///   - construct one canonical representation of Element amongst all equivalent representations
  ///     (i.e. sort sites of a cluster), used for comparison
  /// - std::unique_ptr<SymOpRepresentation> Derived::canonical_transform_impl(Element const  &obj) const
  ///   - access transform that took cluster from unprepared to prepared state
  /// - bool Derived::compare_impl(const Element &A, const Element &B) const
  ///   - used to identify unique equivalents
  /// - InvariantsType Derived::make_invariants_impl(const Element &obj) const
  ///   - construct invariants object for this type of symmetry
  /// - bool Derived::invariants_compare_impl(const Element &A, const Element &B) const
  ///   - first order comparison of orbits
  ///   - speeds up identifying orbits that might contain an element
  /// - void Derived::copy_apply_impl(SymOp const &op, Element obj) const
  ///   - defines application of SymOp to object of type Element
  ///   - allows flexibility for working outside the typical MasterSymGroup/SymGroupRep framework
  ///     or when symmetry application is ambiguous (like applying symmetry to a Eigen::VectorXd)
  ///
  /// Derived may optionally implement:
  /// - bool Derived::inter_orbit_compare_impl const
  ///   - default uses this->invariants_compare and this->compare
  ///
  template<typename Base>
  class SymCompare : public Base {

  public:

    using MostDerived = typename Base::MostDerived;
    using Base::derived;

    using Element = typename traits<MostDerived>::Element;

    /// \brief Applies SymOp to Element
    Element copy_apply(SymOp const &op, Element obj) const;

    /// \brief Prepare an element for comparison via an isometric affine transformation
    Element spatial_prepare(Element obj) const;

    /// \brief Prepare an element for comparison via transformation of its internal representation
    Element representation_prepare(Element obj) const;

    /// \brief Prepare an element for comparison via representation_prepare(), followed by spatial_prepare()
    Element prepare(Element obj) const;

    /// \brief Orders 'prepared' elements
    bool compare(const Element &A, const Element &B) const;

    /// \brief Check equivalence of 'prepared' elements
    bool equal(const Element &A, const Element &B) const;

    // TODO: specify InvariantsType in the SymCompareType instead of traits<Element>
    //
    /// \brief type `InvariantsType` is used to screen elements for equivalence and sort orbits
    using InvariantsType = typename traits<MostDerived>::InvariantsType;

    /// \brief Make orbit invariants from one element in the orbit
    InvariantsType make_invariants(const Element &element) const;

    /// \brief Defines an order for orbits that have different invariants
    bool invariants_compare(const InvariantsType &A, const InvariantsType &B) const;

    /// \brief Check equivalence of invariants
    bool invariants_equal(const InvariantsType &A, const InvariantsType &B) const;

    /// \brief Orders orbit prototypes, breaking invariants_compare ties
    bool inter_orbit_compare(
      const Element &A,
      const InvariantsType &A_invariants,
      const Element &B,
      const InvariantsType &B_invariants) const;

    /// \brief Check equivalence of prototypes in different orbit
    bool inter_orbit_equal(
      const Element &A,
      const InvariantsType &A_invariants,
      const Element &B,
      const InvariantsType &B_invariants) const;

    /// \brief Access transform that took cluster from unprepared to prepared state
    std::unique_ptr<SymOpRepresentation> canonical_transform(Element const  &obj) const;

    /// \brief Access spatial transform that was used during most recent spatial preparation of an element
    SymOp const &spatial_transform() const;

  protected:

    /// \brief Orders orbit prototypes, breaking invariants_compare ties
    bool inter_orbit_compare_impl(
      const Element &A,
      const InvariantsType &A_invariants,
      const Element &B,
      const InvariantsType &B_invariants) const;

  };


  /// \brief Applies SymOp to Element. Default performs standalone copy_apply
  template<typename Base>
  typename SymCompare<Base>::Element SymCompare<Base>::copy_apply(SymOp const &op, Element obj) const {
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
  template<typename Base>
  typename SymCompare<Base>::Element SymCompare<Base>::spatial_prepare(Element obj) const {
    return derived().spatial_prepare_impl(obj);
  }

  /// \brief Prepare an element for comparison via transformation of its internal representation
  ///
  /// - For instance, sort sites of a cluster so comparison may be
  ///   performed more efficiently
  ///
  template<typename Base>
  typename SymCompare<Base>::Element SymCompare<Base>::representation_prepare(Element obj) const {
    return derived().representation_prepare_impl(obj);
  }

  /// \brief Prepare an element for comparison via representation_prepare(), followed by spatial_prepare()
  ///
  template<typename Base>
  typename SymCompare<Base>::Element SymCompare<Base>::prepare(Element obj) const {
    return spatial_prepare(representation_prepare(obj));
  }

  /// \brief Orders 'prepared' elements
  ///
  /// - Returns 'true' to indicate A < B
  /// - Assumes elements are 'prepared' before being compared
  template<typename Base>
  bool SymCompare<Base>::compare(const Element &A, const Element &B) const {
    return derived().compare_impl(A, B);
  }

  /// \brief Check equivalence of 'prepared' elements
  ///
  /// \returns \code !compare(A,B) && !compare(B,A) \endcode
  ///
  /// - Assumes elements are 'prepared' before being compared
  template<typename Base>
  bool SymCompare<Base>::equal(const Element &A, const Element &B) const {
    return !compare(A, B) && !compare(B, A);
  }

  /// \brief Get orbit invariants from one element in the orbit
  template<typename Base>
  typename SymCompare<Base>::InvariantsType SymCompare<Base>::make_invariants(const Element &obj) const {
    return derived().make_invariants_impl(obj);
  }

  /// \brief Defines an order for elements that have different invariants
  ///
  /// - Returns 'true' to indicate A < B
  /// - Equivalence is indicated by \code !invariants_compare(A,B) && !invariants_compare(B,A) \endcode
  /// - Assumes elements are 'prepared' before being compared
  template<typename Base>
  bool SymCompare<Base>::
  invariants_compare(const InvariantsType &A, const InvariantsType &B) const {
    return derived().invariants_compare_impl(A, B);
  }

  /// \brief Check equivalence of invariants
  ///
  /// \returns \code !invariants_compare(A,B) && !invariants_compare(B,A) \endcode
  template<typename Base>
  bool SymCompare<Base>::invariants_equal(const InvariantsType &A, const InvariantsType &B) const {
    return !invariants_compare(A, B) && !invariants_compare(B, A);
  }

  /// \brief Orders orbit prototypes, breaking invariants_compare ties
  ///
  /// - Returns 'true' to indicate A < B
  /// - Equivalence is indicated by \code !inter_orbit_compare(A,B) && !inter_orbit_compare(B,A) \endcode
  /// - Assumes elements are in canonical form
  template<typename Base>
  bool SymCompare<Base>::inter_orbit_compare(
    const Element &A,
    const InvariantsType &A_invariants,
    const Element &B,
    const InvariantsType &B_invariants) const {
    return derived().inter_orbit_compare_impl(A, A_invariants, B, B_invariants);
  }

  /// \brief Check equivalence of prototypes in different orbit
  ///
  /// \returns \code !inter_orbit_compare(A,B) && !inter_orbit_compare(B,A) \endcode
  ///
  /// - Assumes elements are in canonical form
  template<typename Base>
  bool SymCompare<Base>::inter_orbit_equal(
    const Element &A,
    const InvariantsType &A_invariants,
    const Element &B,
    const InvariantsType &B_invariants) const {
    return !inter_orbit_compare(A, A_invariants, B, B_invariants) &&
           !inter_orbit_compare(B, B_invariants, A, A_invariants);
  }

  /// \brief Access transform that took cluster from unprepared to prepared state
  template<typename Base>
  std::unique_ptr<SymOpRepresentation> SymCompare<Base>::canonical_transform(Element const  &obj) const {
    return derived().canonical_transform_impl(obj);
  }

  /// \brief Access spatial transform that was used during most recent spatial preparation of an element
  template<typename Base>
  SymOp const &SymCompare<Base>::spatial_transform() const {
    return derived().spatial_transform_impl();
  }

  /// \brief Orders orbit prototypes, breaking invariants_compare ties
  ///
  /// - Returns 'true' to indicate A < B
  /// - Equivalence is indicated by \code !compare(A,B) && !compare(B,A) \endcode
  ///
  /// Implementation:
  /// - First, check invariants_compare
  /// - Break ties with, compare(A, B)
  template<typename Base>
  bool SymCompare<Base>::inter_orbit_compare_impl(
    const Element &A,
    const InvariantsType &A_invariants,
    const Element &B,
    const InvariantsType &B_invariants) const {

    // first compare invariants
    if(this->invariants_compare(A_invariants, B_invariants)) {
      return true;
    }
    if(this->invariants_compare(B_invariants, A_invariants)) {
      return false;
    }

    // next compare A and B
    return this->compare(A, B);
  }

}

#endif
