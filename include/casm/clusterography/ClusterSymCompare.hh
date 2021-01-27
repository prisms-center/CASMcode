#ifndef CASM_ClusterSymCompare
#define CASM_ClusterSymCompare

#include <memory>

#include "casm/clusterography/ClusterSymCompareDecl.hh"
#include "casm/crystallography/IntegralCoordinateWithin.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/global/eigen.hh"
#include "casm/symmetry/SymCompare.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/SymTools.hh"

namespace CASM {

class Structure;

/** \defgroup ClusterSymCompare

    \brief Functions and classes for creating cluster SymCompare functors for
   orbit generation \ingroup Clusterography \ingroup IntegralCluster

*/

/// \brief CRTP Base class for common cluster symmetry and comparison
/// implementations
///
/// ClusterSymCompare can be used for common implementations of SymCompareType
/// methods for cluster orbits. To do so, it requires that a traits class
/// (traits<SymCompareType>) is implemented with:
/// - typedef Element; // type of Cluster, ex: IntegralCluster
/// - typedef InvariantsType; // ex: ClusterInvariants
/// - static Element copy_apply(SymOp const&, Element, SymCompareType const&);
/// - static InvariantsType make_invariants(Element const&, SymCompareType
/// const&);
///
/// The orbit element type, Element, is required to implement:
/// - 'Element::operator<(Element B)'
/// - 'Permutation Element::sort_permutation()'
///
/// The orbit invariants type, InvariantsType, is required to implement:
/// - 'CASM::compare(InvariantsType, InvariantsType, double tol)'
///
/// ClusterSymCompare itself implements:
/// - 'make_invariants_impl':
///   - uses 'traits<SymCompareType>::make_invariants(Element, SymCompareType)'
/// - 'invariants_compare_impl'
///   - uses 'CASM::compare(InvariantsType, InvariantsType, double tol)';
/// - 'compare_impl':
///   - uses 'Element::operator<(Element B)'
/// - 'copy_apply_impl':
///   - implemented by SymCompareType
///   - uses 'traits<SymCompareType>::copy_apply(SymOp, Element,
///   SymCompareType)'
/// - 'canonical_transform_impl':
///   - uses 'Element::sort_permutation'
/// - 'inter_orbit_compare_impl':
///   - uses the default implemented by SymCompare
///   - uses invariants_compare_impl and compare_impl
///
/// The methods ClusterSymCompare implements are overriden if an implementation
/// exists in final SymCompareType implementations. For example, this is done by
/// WithinScelSymCompare to ensure all sites are within the supercell before
/// finding the canonical transform permutation.
///
/// Final SymCompareType implementations inheriting from ClusterSymCompare still
/// require:
/// - spatial_prepare_impl
///   - typically uses traits<SymCompareType>::position
/// - spatial_transform_impl
/// - representation_prepare_impl
///   - typically uses Element::sort
///
/// \ingroup Clusterography
///
template <typename Base>
class ClusterSymCompare : public Base {
 public:
  typedef typename Base::MostDerived MostDerived;
  using typename Base::Element;  /// Element (of orbit) refers to Cluster, not
                                 /// element of Cluster
  typedef Element ClusterType;
  using Base::derived;
  using typename Base::InvariantsType;

 protected:
  /// \brief Make orbit invariants from one element in the orbit
  InvariantsType make_invariants_impl(const ClusterType &obj) const;

  /// \brief Orders 'prepared' elements in the same orbit
  bool invariants_compare_impl(const InvariantsType &A,
                               const InvariantsType &B) const;

  /// \brief Compares 'prepared' clusters
  bool compare_impl(const ClusterType &A, const ClusterType &B) const;

  /// \brief Applies SymOp to cluster
  ClusterType copy_apply_impl(SymOp const &op, ClusterType obj) const;

  /// \brief Returns transformation that takes 'obj' to its prepared (canonical)
  /// form
  // For now, this is the the sorting permutation
  std::unique_ptr<SymOpRepresentation> canonical_transform_impl(
      ClusterType const &obj) const;
};

/* -- AperiodicSymCompare<IntegralCluster> Declaration
 * ------------------------------------- */

/// \brief Comparisons of clusters using aperiodic symmetry
///
/// Before doing a comparison, AperiodicSymCompare does not spatially transform
/// the cluster in any way, as opposed to PrimPeriodSymCompare which translates
/// the cluster so that the first site is within the origin unit cell.
///
/// Beyond ClusterSymCompare requirements, requires:
/// - 'Element& Element::sort();'
///
/// \ingroup Clusterography
/// \ingroup ClusterSymCompare
///
template <typename Element>
class AperiodicSymCompare<Element>
    : public ClusterSymCompare<
          SymCompare<CRTPBase<AperiodicSymCompare<Element>>>> {
 public:
  typedef Structure PrimType;
  typedef std::shared_ptr<PrimType const> PrimType_ptr;

  /// \brief Constructor
  ///
  /// \param prim_ptr Prim structure
  /// \param tol Tolerance for invariants_compare of site-to-site distances
  ///
  AperiodicSymCompare(PrimType_ptr prim_ptr, double tol);

  const PrimType &prim() const { return *m_prim; }

  /// \brief Return tolerance
  double tol() const { return this->m_tol; }

 private:
  typedef ClusterSymCompare<SymCompare<CRTPBase<AperiodicSymCompare<Element>>>>
      Base;
  friend traits<AperiodicSymCompare<Element>>;
  friend Base;
  friend SymCompare<CRTPBase<AperiodicSymCompare<Element>>>;

  /// \brief Prepare an element for comparison via an isometric affine
  /// transformation
  ///
  /// - For aperiodic cases, no isometric transformations are allowed, so apply
  /// and return identity
  Element spatial_prepare_impl(Element obj) const;

  /// \brief Access spatial transform that was used during most recent spatial
  /// preparation of an element
  /// - Always identity
  SymOp const &spatial_transform_impl() const;

  /// \brief Prepare an element for comparison via transformation of its
  /// internal representation
  ///
  /// - Returns sorted
  Element representation_prepare_impl(Element obj) const;

  PrimType_ptr m_prim;

  double m_tol;

  /// \brief Spatial transform that reproduces most recent application of
  /// SymCompare::spatial_prepare()
  /// - Default SymOp constructor initializes to identity
  mutable SymOp m_spatial_transform;
};

/// \brief Comparisons of clusters using prim periodic symmetry
///
/// Before doing a comparison of cluster sites, PrimPeriodicSymCompare
/// translates the cluster so that the first site is within the origin unit
/// cell.
///
/// Beyond ClusterSymCompare requirements, requires:
/// - Element must be translatable (Element + UnitCell, Element - Unitcell are
/// valid)
/// - 'static xtal::UnitCellCoord
/// traits<PrimPeriodicSymCompare<Element>>::position(
///     Element const&,
///     PrimPeriodicSymCompare<Element> const& sym_compare);
/// - 'Element& Element::sort();'
///
/// \ingroup Clusterography
/// \ingroup ClusterSymCompare
///
template <typename Element>
class PrimPeriodicSymCompare<Element>
    : public ClusterSymCompare<
          SymCompare<CRTPBase<PrimPeriodicSymCompare<Element>>>> {
 public:
  typedef Structure PrimType;
  typedef std::shared_ptr<PrimType const> PrimType_ptr;

  /// \brief Constructor
  ///
  /// \param prim_ptr Prim structure
  /// \param tol Tolerance for invariants_compare of site-to-site distances
  ///
  PrimPeriodicSymCompare(PrimType_ptr prim_ptr, double tol);

  const PrimType &prim() const { return *m_prim; }

  /// \brief Return tolerance
  double tol() const { return this->m_tol; }

 private:
  typedef ClusterSymCompare<
      SymCompare<CRTPBase<PrimPeriodicSymCompare<Element>>>>
      Base;
  friend traits<PrimPeriodicSymCompare<Element>>;
  friend Base;
  friend SymCompare<CRTPBase<PrimPeriodicSymCompare<Element>>>;

  /// \brief Prepare an element for comparison via an isometric affine
  /// transformation
  ///
  /// - Applies lattice translation such that first site of cluster is in
  /// UnitCell (0,0,0)
  Element spatial_prepare_impl(Element obj) const;

  /// \brief Access spatial transform that was used during most recent spatial
  /// preparation of an element
  SymOp const &spatial_transform_impl() const;

  /// \brief Prepare an element for comparison via transformation of its
  /// internal representation
  ///
  /// - Returns sorted
  Element representation_prepare_impl(Element obj) const;

  /// Pointer to the primitive structure, necessary to apply symmetry to the
  /// Element
  PrimType_ptr m_prim;

  double m_tol;

  /// \brief Spatial transform that reproduces most recent application of
  /// SymCompare::spatial_prepare()
  /// - Default SymOp constructor initializes to identity
  mutable SymOp m_spatial_transform;
};

/// \brief Comparisons of clusters using supercell periodic symmetry, but
/// without periodic images
///
/// Before doing a comparison, ScelPeriodicSymCompare translates the cluster so
/// that the first site is within the supercell, as opposed to
/// WithinScelSymCompare which moves all sites in the cluster within the
/// supercell.
///
/// ScelPeriodicSymCompare uses the direct distance between sites for orbit
/// invariants without accounting for periodic images, while
/// WithinScelSymCompare uses supercell periodic image minimum distance for the
/// orbit invariants.
///
/// Beyond ClusterSymCompare requirements, requires:
/// - Element must be translatable (Element + UnitCell, Element - Unitcell are
/// valid)
/// - static xtal::UnitCellCoord
/// traits<PrimPeriodicSymCompare<Element>>::position(
///     Element const&,
///     PrimPeriodicSymCompare<Element> const& sym_compare);
/// - Element& Element::sort();
///
/// \ingroup Clusterography
/// \ingroup ClusterSymCompare
///
template <typename Element>
class ScelPeriodicSymCompare<Element>
    : public ClusterSymCompare<
          SymCompare<CRTPBase<ScelPeriodicSymCompare<Element>>>> {
 public:
  typedef Structure PrimType;
  typedef std::shared_ptr<PrimType const> PrimType_ptr;
  typedef Eigen::Matrix3l transf_mat_type;

  /// \brief Constructor
  ///
  /// \param prim_ptr Prim structure
  /// \param transf_mat Prim to supercell transformation matrix
  /// \param tol Tolerance for invariants_compare of site-to-site distances
  ///
  ScelPeriodicSymCompare(PrimType_ptr prim_ptr, transf_mat_type transf_mat,
                         double tol);

  const PrimType &prim() const { return *m_prim; }

  /// Prim to supercell transformation matrix
  transf_mat_type const &transf_mat() const { return m_transf_mat; }

  /// Return tolerance
  double tol() const { return this->m_tol; }

 private:
  typedef ClusterSymCompare<
      SymCompare<CRTPBase<ScelPeriodicSymCompare<Element>>>>
      Base;
  friend traits<ScelPeriodicSymCompare<Element>>;
  friend Base;
  friend SymCompare<CRTPBase<ScelPeriodicSymCompare<Element>>>;

  /// \brief Prepare an element for comparison via an isometric affine
  /// transformation
  ///
  /// - Applies superlattice translation such that first site of cluster is
  /// within supercell
  Element spatial_prepare_impl(Element obj) const;

  /// \brief Access spatial transform that was used during most recent spatial
  /// preparation of an element
  SymOp const &spatial_transform_impl() const;

  /// \brief Prepare an element for comparison via transformation of its
  /// internal representation
  ///
  /// - Returns sorted
  Element representation_prepare_impl(Element obj) const;

  /// Pointer to the primitive structure, necessary to apply symmetry to the
  /// Element
  PrimType_ptr m_prim;

  /// Prim to supercell transformation matrix
  transf_mat_type m_transf_mat;

  /// Bring UnitCellCoord within the supercell
  xtal::IntegralCoordinateWithin_f m_bring_within_f;

  double m_tol;

  /// \brief Spatial transform that reproduces most recent application of
  /// SymCompare::spatial_prepare()
  /// - Default SymOp constructor initializes to identity
  mutable SymOp m_spatial_transform;
};

// /// \brief Comparisons of clusters using supercell periodic symmetry, with
// periodic images
// ///
// /// Before doing a comparison, WithinScelSymCompare moves all sites in the
// cluster within the
// /// supercell, as opposed to ScelPeriodicSymCompare which only translates the
// cluster so that the
// /// first site is within the supercell.
// ///
// /// WithinScelSymCompare uses supercell periodic image minimum distance for
// the orbit
// /// invariants, while ScelPeriodicSymCompare just uses the direct distance
// without accounting for
// /// periodic images.
// ///
// /// To implement, traits<WithinScelSymCompare<Element>> is required to have:
// /// - static IntegralCluster bring_within(
// ///     IntegralCluster clust,
// ///     WithinScelSymCompare<IntegralCluster> const& sym_compare);
// /// - Element& Element::sort();
// ///
// /// \ingroup Clusterography
// /// \ingroup ClusterSymCompare
// ///
// template <typename Element>
// class WithinScelSymCompare<Element>:
//   public
//   ClusterSymCompare<SymCompare<CRTPBase<WithinScelSymCompare<Element>>>> {
//
// public:
//   typedef Structure PrimType;
//   typedef std::shared_ptr<PrimType const> PrimType_ptr;
//   typedef Eigen::Matrix3l transf_mat_type;
//
//   /// \brief Constructor
//   ///
//   /// \param prim_ptr Prim structure
//   /// \param transf_mat Prim to supercell transformation matrix
//   /// \param tol Tolerance for invariants_compare of site-to-site distances
//   ///
//   WithinScelSymCompare(
//     PrimType_ptr prim_ptr,
//     transf_mat_type transf_mat,
//     double tol);
//
//   const PrimType &prim() const {
//     return *m_prim;
//   }
//
//   /// Prim to supercell transformation matrix
//   transf_mat_type const &transf_mat() const {
//     return m_transf_mat;
//   }
//
//   /// \brief Return tolerance
//   double tol() const {
//     return this->m_tol;
//   }
//
// private:
//   typedef
//   ClusterSymCompare<SymCompare<CRTPBase<WithinScelSymCompare<Element>>>>
//   Base; friend traits<WithinScelSymCompare<Element>>; friend
//   SymCompare<CRTPBase<WithinScelSymCompare<Element>>>; friend Base;
//
//   /// \brief Returns transformation that takes 'obj' to its prepared
//   (canonical) form
//   ///
//   std::unique_ptr<SymOpRepresentation> canonical_transform_impl(Element const
//   &obj) const;
//
//   /// \brief Prepare an element for comparison via an isometric affine
//   transformation
//   ///
//   /// - Applies superlattice translation such that first site of cluster is
//   within supercell Element spatial_prepare_impl(Element obj) const;
//
//   /// \brief Access spatial transform that was used during most recent
//   spatial preparation of an element
//   /// - Always identity
//   SymOp const &spatial_transform_impl() const;
//
//   /// \brief Prepare an element for comparison via transformation of its
//   internal representation
//   ///
//   /// - Returns sorted
//   Element representation_prepare_impl(Element obj) const;
//
//   /// Pointer to the primitive structure, necessary to apply symmetry to the
//   Element PrimType_ptr m_prim;
//
//   /// Prim to supercell transformation matrix
//   transf_mat_type m_transf_mat;
//
//   /// Bring UnitCellCoord within the supercell
//   xtal::IntegralCoordinateWithin_f m_bring_within_f;
//
//   double m_tol;
//
//   /// \brief Spatial transform that reproduces most recent application of
//   SymCompare::spatial_prepare()
//   /// - Default SymOp constructor initializes to identity
//   mutable SymOp m_spatial_transform;
// };

}  // namespace CASM

#endif
