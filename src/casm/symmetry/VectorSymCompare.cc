#include "casm/symmetry/VectorSymCompare_impl.hh"

namespace CASM {

  VectorInvariants::VectorInvariants(Eigen::VectorXcd const &vector):
    m_cols(vector.cols()), m_norm(vector.norm()) {}

  double VectorInvariants::cols() const {
    return m_cols;
  }

  double VectorInvariants::norm() const {
    return m_norm;
  }

  bool almost_equal(VectorInvariants const &A_invariants, VectorInvariants const &B_invariants, double tol) {
    return A_invariants.cols() == B_invariants.cols()
           && almost_equal(A_invariants.norm(), B_invariants.norm(), tol);
  }

  bool compare(VectorInvariants const &A_invariants, VectorInvariants const &B_invariants, double tol) {
    if(A_invariants.cols() == B_invariants.cols()) {
      return CASM::compare(A_invariants.norm(), B_invariants.norm(), tol);
    }
    return A_invariants.cols() < B_invariants.cols();
  }

}
