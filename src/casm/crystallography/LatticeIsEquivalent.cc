#include "casm/crystallography/LatticeIsEquivalent.hh"

#include "casm/symmetry/SymOp.hh"

namespace CASM {

  LatticeIsEquivalent::LatticeIsEquivalent(const Lattice &lat, double _tol) :
    m_lat(lat), m_tol(_tol) {}


  /// Is this lattice the same, even if they have different lattice vectors
  bool LatticeIsEquivalent::operator()(const Lattice &B) const {
    Eigen::Matrix3d T = lat_column_mat().inverse() * B.lat_column_mat();
    return is_unimodular(T, m_tol);
  }

  /// Is this lattice equivalent to apply(op, *this)
  bool LatticeIsEquivalent::operator()(const SymOp &op) const {
    return (*this)(op.matrix());
  }

  /// Is this lattice equivalent to apply(op, *this)
  bool LatticeIsEquivalent::operator()(const Eigen::Matrix3d &cart_op) const {
    Eigen::Matrix3d tfrac_op;

    tfrac_op = lat_column_mat().inverse() * cart_op * lat_column_mat();

    //Use a soft tolerance of 1% to see if further screening should be performed
    if(!almost_equal(1.0, std::abs(tfrac_op.determinant()), 0.01) || !is_integer(tfrac_op, 0.01)) {
      return false;
    }

    //make tfrac_op integer.
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        tfrac_op(i, j) = round(tfrac_op(i, j));
      }
    }

    return _check(tfrac_op);
  }

  /// Is this lattice equivalent to apply(op, *this)
  bool LatticeIsEquivalent::operator()(const Eigen::Matrix3i &tfrac_op) const {

    //false if determinant is not 1, because it doesn't preserve volume
    if(std::abs(tfrac_op.determinant()) != 1) {
      return false;
    }

    return _check(tfrac_op.cast<double>());
  }

  const Lattice &LatticeIsEquivalent::lat() const {
    return m_lat;
  }

  double LatticeIsEquivalent::map_error() const {
    return m_map_error;
  }

  Eigen::Matrix3d LatticeIsEquivalent::cart_op() const {
    return m_cart_op;
  }

  SymOp LatticeIsEquivalent::sym_op() const {
    return SymOp(cart_op(), map_error());
  }

  ///Find the effect of applying symmetry to the lattice vectors
  bool LatticeIsEquivalent::_check(const Eigen::Matrix3d &tfrac_op) const {

    // If symmetry is perfect, then ->  cart_op * lat_column_mat() == lat_column_mat() * frac_op  by definition
    // If we assum symmetry is imperfect, then ->   cart_op * lat_column_mat() == F * lat_column_mat() * frac_op
    // where 'F' is the displacement gradient tensor imposed by frac_op
    m_cart_op = lat_column_mat() * tfrac_op * inv_lat_column_mat();

    // tMat uses some matrix math to get F.transpose()*F*lat_column_mat();
    Eigen::Matrix3d tMat = m_cart_op.transpose() * lat_column_mat() * tfrac_op;

    // Subtract lat_column_mat() from tMat, leaving us with (F.transpose()*F - Identity)*lat_column_mat().
    // This is 2*E*lat_column_mat(), where E is the green-lagrange strain
    tMat = (tMat - lat_column_mat()) / 2.0;

    //... and then multiplying by the transpose...
    tMat = tMat * tMat.transpose();

    // The diagonal elements of tMat describe the square of the distance by which the transformed vectors 'miss' the original vectors
    if(tMat(0, 0) < m_tol * m_tol && tMat(1, 1) < m_tol * m_tol && tMat(2, 2) < m_tol * m_tol) {
      m_map_error = sqrt(tMat.diagonal().maxCoeff());
      return true;
    }
    return false;
  }

  const Eigen::Matrix3d &LatticeIsEquivalent::lat_column_mat() const {
    return m_lat.lat_column_mat();
  }

  const Eigen::Matrix3d &LatticeIsEquivalent::inv_lat_column_mat() const {
    return m_lat.inv_lat_column_mat();
  }

}

