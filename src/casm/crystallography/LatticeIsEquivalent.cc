#include "casm/crystallography/LatticeIsEquivalent.hh"

#include "casm/crystallography/SymTools.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
namespace xtal {

LatticeIsEquivalent::LatticeIsEquivalent(const Lattice &_lat) : m_lat(_lat) {}

/// Checks if lat = other*U, with unimodular U
bool LatticeIsEquivalent::operator()(const Lattice &other, double tol) const {
  tol = tol < 0 ? m_lat.tol() : tol;
  m_U = other.lat_column_mat().inverse() * m_lat.lat_column_mat();
  return is_unimodular(m_U, tol);
}

/// Returns U found for last check
Eigen::Matrix3d LatticeIsEquivalent::U() const { return m_U; }

IsPointGroupOp::IsPointGroupOp(const Lattice &lat) : m_lat(lat) {}

/// Is this lattice equivalent to apply(op, *this)
bool IsPointGroupOp::operator()(const SymOp &op) const {
  return (*this)(get_matrix(op));
}

/// Is this lattice equivalent to apply(op, *this)
bool IsPointGroupOp::operator()(const Eigen::Matrix3d &cart_op) const {
  Eigen::Matrix3d tfrac_op;

  tfrac_op = lat_column_mat().inverse() * cart_op * lat_column_mat();

  // Use a soft tolerance of 1% to see if further screening should be performed
  if (!almost_equal(1.0, std::abs(tfrac_op.determinant()), 0.01) ||
      !is_integer(tfrac_op, 0.01)) {
    return false;
  }

  return _check(round(tfrac_op));
}

/// Is this lattice equivalent to apply(op, *this)
bool IsPointGroupOp::operator()(const Eigen::Matrix3i &tfrac_op) const {
  // false if determinant is not 1, because it doesn't preserve volume
  if (std::abs(tfrac_op.determinant()) != 1) {
    return false;
  }

  return _check(tfrac_op.cast<double>());
}

double IsPointGroupOp::map_error() const { return m_map_error; }

Eigen::Matrix3d IsPointGroupOp::cart_op() const { return m_cart_op; }

SymOp IsPointGroupOp::sym_op() const {
  return SymOp::point_operation(cart_op());
}

/// Find the effect of applying symmetry to the lattice vectors
bool IsPointGroupOp::_check(const Eigen::Matrix3d &tfrac_op) const {
  // If symmetry is perfect, then ->  cart_op * lat_column_mat() ==
  // lat_column_mat() * frac_op  by definition If we assum symmetry is
  // imperfect, then ->   cart_op * lat_column_mat() == F * lat_column_mat() *
  // frac_op where 'F' is the displacement gradient tensor imposed by frac_op
  m_cart_op = lat_column_mat() * tfrac_op * inv_lat_column_mat();

  // tMat uses some matrix math to get F.transpose()*F*lat_column_mat();
  Eigen::Matrix3d tMat = m_cart_op.transpose() * lat_column_mat() * tfrac_op;

  // Subtract lat_column_mat() from tMat, leaving us with (F.transpose()*F -
  // Identity)*lat_column_mat(). This is 2*E*lat_column_mat(), where E is the
  // green-lagrange strain
  tMat = (tMat - lat_column_mat()) / 2.0;

  //... and then multiplying by the transpose...
  tMat = tMat * tMat.transpose();

  // The diagonal elements of tMat describe the square of the distance by which
  // the transformed vectors 'miss' the original vectors

  double sqr_tol = m_lat.tol() * m_lat.tol();
  if (tMat(0, 0) < sqr_tol && tMat(1, 1) < sqr_tol && tMat(2, 2) < sqr_tol) {
    m_map_error = sqrt(tMat.diagonal().maxCoeff());
    return true;
  }
  return false;
}

const Eigen::Matrix3d &IsPointGroupOp::lat_column_mat() const {
  return m_lat.lat_column_mat();
}

const Eigen::Matrix3d &IsPointGroupOp::inv_lat_column_mat() const {
  return m_lat.inv_lat_column_mat();
}

bool is_equivalent(const Lattice &ref_lattice, const Lattice &other) {
  LatticeIsEquivalent is_equiv(ref_lattice);
  return is_equiv(other);
}

}  // namespace xtal
}  // namespace CASM
