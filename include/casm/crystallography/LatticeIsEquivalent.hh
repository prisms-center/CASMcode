#ifndef CASM_LatticeIsEquivalent
#define CASM_LatticeIsEquivalent

#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/external/Eigen/Dense"

namespace CASM {
namespace xtal {

struct SymOp;

/// \brief Lattice comparisons
///
/// Does comparisons of the form:
///     copy_apply(A, lat) ?= copy_apply(B, other)*U,
///
///   where lat and other are lattices, represented as column matrices
///   A & B are symmetry operations
///   U is a unimodular (3x3, integer, w/ abs(det(T))==1) transformation matrix
///
/// \ingroup Lattice
/// \ingroup IsEquivalent
///
class LatticeIsEquivalent {
 public:
  LatticeIsEquivalent(const Lattice &_lat);

  /// Checks if lat = other*U, with unimodular U, tolerance is relative,
  /// default tolerance is sentinal value to default to  m_lat.tol()
  bool operator()(const Lattice &other, double tol = -1.) const;

  /// Returns U found for last check
  Eigen::Matrix3d U() const;

 private:
  Lattice m_lat;
  mutable Eigen::Matrix3d m_U;
};

/// Checks if operations are point group operations
class IsPointGroupOp {
 public:
  IsPointGroupOp(const Lattice &lat);

  /// Checks if ref_lat = cart_op*ref_lat*transf_mat(), for any transf_mat()
  bool operator()(const SymOp &cart_op) const;

  template <typename ExternSymOp>
  bool operator()(const ExternSymOp &cart_op) const {
    return this->operator()(adapter::Adapter<SymOp, ExternSymOp>()(cart_op));
  }

  /// Checks if ref_lat = cart_op*ref_lat*transf_mat(), for any transf_mat()
  bool operator()(const Eigen::Matrix3d &cart_op) const;

  /// Checks if ref_lat = (ref_lat*frac_op)*transf_mat(), for any transf_mat()
  bool operator()(const Eigen::Matrix3i &frac_op) const;

  /// Return the mapping error, calculated after performing an equivalence check
  double map_error() const;

  /// If evaluates true, then ref_lat == cart_op()*L*transf_mat() to the
  /// specified tolerance
  Eigen::Matrix3d cart_op() const;

  SymOp sym_op() const;

 private:
  /// Find the effect of applying symmetry to the lattice vectors
  bool _check(const Eigen::Matrix3d &tfrac_op) const;

  const Eigen::Matrix3d &lat_column_mat() const;

  const Eigen::Matrix3d &inv_lat_column_mat() const;

  Lattice m_lat;
  mutable double m_map_error;
  mutable Eigen::Matrix3d m_cart_op;
};

/// Check if ref_lattice = other*U, where U is unimodular
bool is_equivalent(const Lattice &ref_lattice, const Lattice &other);

}  // namespace xtal
}  // namespace CASM

#endif
