#ifndef CASM_LatticeIsEquivalent
#define CASM_LatticeIsEquivalent

#include "casm/external/Eigen/Dense"
#include "casm/crystallography/Lattice.hh"

namespace CASM {

  /// \brief Putting all the Lattice comparisons in one place
  class LatticeIsEquivalent {

  public:

    LatticeIsEquivalent(const Lattice &lat, double _tol = TOL);


    /// Is this lattice the same, even if they have different lattice vectors
    bool operator()(const Lattice &B) const;

    /// Is this lattice equivalent to apply(op, *this)
    bool operator()(const SymOp &op) const;

    /// Is this lattice equivalent to apply(op, *this)
    bool operator()(const Eigen::Matrix3d &cart_op) const;

    /// Is this lattice equivalent to apply(op, *this)
    bool operator()(const Eigen::Matrix3i &tfrac_op) const;

    const Lattice &lat() const;

    double map_error() const;

    Eigen::Matrix3d cart_op() const;

    SymOp sym_op() const;

  private:

    ///Find the effect of applying symmetry to the lattice vectors
    bool _check(const Eigen::Matrix3d &tfrac_op) const;

    const Eigen::Matrix3d &lat_column_mat() const;

    const Eigen::Matrix3d &inv_lat_column_mat() const;

    Lattice m_lat;
    double m_tol;
    mutable double m_map_error;
    mutable Eigen::Matrix3d m_cart_op;

  };

}

#endif

