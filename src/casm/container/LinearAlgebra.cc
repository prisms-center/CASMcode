#include "casm/container/LinearAlgebra.hh"
#include "casm/container/Counter.hh"

#include "casm/external/boost.hh"

namespace CASM {

  //*******************************************************************************************
  /**
   * Given a matrix, original_mat, this function calculates
   * both the Hermitian and anti-Hermitian parts.
   *
   * @param[in] original_mat The matrix we want to get the
   *                         Hermitian and anti-Hermitian
   *                         parts of
   * @param[in,out] hermitian_mat The Hermitian matrix
   * @param[in,out] antihermitian_mat The anti-Hermitian matrix
   */
  //*******************************************************************************************
  //Added by Ivy 11/28/12
  void get_Hermitian(Eigen::MatrixXcd &original_mat, Eigen::MatrixXcd &hermitian_mat, Eigen::MatrixXcd &antihermitian_mat) {

    Eigen::MatrixXcd conj_trans = original_mat.conjugate().transpose();

    hermitian_mat = 0.5 * (original_mat + conj_trans);
    antihermitian_mat = 0.5 * (original_mat - conj_trans);
  }

  //*******************************************************************************************
  /**
   * Returns true if the specified matrix is Hermitian.
   *
   * Checks if the matrix is the same as the conjugate
   */
  //*******************************************************************************************
  //Added by Ivy 11/28/12
  bool is_Hermitian(Eigen::MatrixXcd &mat) {

    if((mat - mat.conjugate().transpose()).isZero(TOL)) {
      return true;
    }
    return false;
  }



  /// \brief Return the hermite normal form, M == H*V
  ///
  /// \param M a square integer matrix
  ///
  /// \throws std::runtime_error if M is not full rank
  /// \returns std::pair<H,V>
  ///
  /// H set to upper triangular matrix H with:
  /// - H(i,i) > 0
  /// - 0 <= H(r,i) < H(i,i) for r > i
  /// - H(r,i) == 0 for r < i
  /// V set to unimodular matrix V
  /// - V.determinant() = +/- 1
  ///
  std::pair<Eigen::MatrixXi, Eigen::MatrixXi> hermite_normal_form(const Eigen::MatrixXi &M) {

    if(M.rows() != M.cols()) {
      throw std::runtime_error(
        std::string("Error in hermite_normal_form(const Eigen::MatrixXi &M)\n") +
        "  M must be square.");
    }

    if(boost::math::iround(M.cast<double>().determinant()) == 0) {
      throw std::runtime_error(
        std::string("Error in hermite_normal_form(const Eigen::MatrixXi &M)\n") +
        "  M must be full rank.");
    }

    Eigen::MatrixXi V = Eigen::MatrixXi::Identity(M.rows(), M.cols());
    Eigen::MatrixXi H = M;



    int N = M.rows();
    int i, r, c;


    // get non-zero entries along H diagonal
    for(i = N - 1; i >= 0; i--) {

      // swap columns if necessary to get non-zero entry at H(i,i)
      if(H(i, i) == 0) {
        for(c = i - 1; c >= 0; c--) {
          if(H(i, c) != 0) {
            H.col(c).swap(H.col(i));
            V.row(c).swap(V.row(i));
            break;
          }
        }
      }

      // get positive entry at H(i,i)
      if(H(i, i) < 0) {
        H.col(i) *= -1;
        V.row(i) *= -1;
      }

      // do column operations to get zeros to left of H(i,i)
      for(c = i - 1; c >= 0; c--) {
        if(H(i, c) == 0)
          continue;

        while(H(i, c) != 0) {

          while(H(i, c) < 0) {
            H.col(c) += H.col(i);
            V.row(i) -= V.row(c);
          }

          while(H(i, i) <= H(i, c) && H(i, c) != 0) {
            H.col(c) -= H.col(i);
            V.row(i) += V.row(c);
          }

          while(H(i, c) < H(i, i) && H(i, c) != 0) {
            H.col(i) -= H.col(c);
            V.row(c) += V.row(i);
          }

        }
      }
    }

    // now we have M = H*V, where H is upper triangular
    // so use column operations to enforce 0 <= H(i,c) < H(i,i) for c > i

    // columns left to right
    for(c = 1; c < H.cols(); c++) {
      // rows above the diagonal to the top
      for(r = c - 1; r >= 0; r--) {
        while(H(r, c) < 0) {
          H.col(c) += H.col(r);
          V.row(r) -= V.row(c);
        }
        while(H(r, c) >= H(r, r)) {
          H.col(c) -= H.col(r);
          V.row(r) += V.row(c);
        }
      }
    }

    // now we have M = H*V
    return std::make_pair(H, V);
  }


  /// \brief Get angle, in radians, between two vectors on range [0,pi]
  double angle(const Eigen::Ref<const Eigen::Vector3d> &a, const Eigen::Ref<const Eigen::Vector3d> &b) {
    return acos(a.dot(b) / (a.norm() * b.norm()));
  }

  ///return signed angle, in radians, between -pi and pi that describe separation in direction of two vectors
  double signed_angle(const Eigen::Ref<const Eigen::Vector3d> &a,
                      const Eigen::Ref<const Eigen::Vector3d> &b,
                      const Eigen::Ref<const Eigen::Vector3d> &pos_ref) {
    if(pos_ref.dot(a.cross(b)) < 0) {
      return -angle(a, b);
    }
    else
      return angle(a, b);
  }

  /// \brief Round entries that are within tol of being integer to that integer value
  Eigen::MatrixXd pretty(const Eigen::MatrixXd &M, double tol) {
    Eigen::MatrixXd Mp(M);
    for(int i = 0; i < M.rows(); i++) {
      for(int j = 0; j < M.cols(); j++) {
        if(std::abs(std::round(M(i, j)) - M(i, j)) < tol) {
          Mp(i, j) = std::round(M(i, j));
        }
      }
    }
    return Mp;
  }

  std::vector<Eigen::Matrix3i> _unimodular_matrices(bool positive, bool negative) {
    std::vector<Eigen::Matrix3i> uni_det_mats;
    int totalmats = 3480;

    if(positive && negative) {
      totalmats = totalmats * 2;
    }

    uni_det_mats.reserve(totalmats);

    EigenCounter<Eigen::Matrix3i> transmat_count(Eigen::Matrix3i::Constant(-1), Eigen::Matrix3i::Constant(1), Eigen::Matrix3i::Constant(1));

    for(; transmat_count.valid(); ++transmat_count) {
      if(positive && transmat_count.current().determinant() == 1) {
        uni_det_mats.push_back(transmat_count.current());
      }

      if(negative && transmat_count.current().determinant() == -1) {
        uni_det_mats.push_back(transmat_count.current());
      }
    }

    return uni_det_mats;
  }

  const std::vector<Eigen::Matrix3i> &positive_unimodular_matrices() {
    static std::vector<Eigen::Matrix3i> static_positive(_unimodular_matrices(true, false));
    return static_positive;
  }

  const std::vector<Eigen::Matrix3i> &negative_unimodular_matrices() {
    static std::vector<Eigen::Matrix3i> static_negative(_unimodular_matrices(true, false));
    return static_negative;
  }

  const std::vector<Eigen::Matrix3i> &unimodular_matrices() {
    static std::vector<Eigen::Matrix3i> static_all(_unimodular_matrices(true, false));
    return static_all;
  }
}

