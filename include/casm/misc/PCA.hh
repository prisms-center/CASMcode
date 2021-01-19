#ifndef CASM_PCA
#define CASM_PCA

#include "casm/external/Eigen/Dense"

namespace CASM {

/// \brief Principle component analysis
class PCA {
 public:
  typedef Eigen::MatrixXd::Index size_type;

  /// \brief Constructor accepting a column vector matrix M containing points
  /// [input...]
  template <typename Derived>
  PCA(const Eigen::MatrixBase<Derived> &M, double singular_value_tol = 1e-14) {
    // set the mean to zero
    Eigen::MatrixXd mat_mean_zero = M;
    mat_mean_zero.colwise() -= (M.rowwise().mean());

    // principle component analysis to rotate to input space that to its range
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat_mean_zero, Eigen::ComputeFullU);
    int rank = 0;
    for (Index i = 0; i < svd.singularValues().size(); i++) {
      if (std::abs(svd.singularValues()[i]) <= singular_value_tol) {
        break;
      } else {
        rank++;
      }
    }

    m_reduce = svd.matrixU().block(0, 0, M.rows(), rank).transpose();
  }

  /// \brief Initial dimension, equivalent to pca.reduce().cols()
  size_type dim() const { return m_reduce.cols(); }

  /// \brief Rank of the input, equivalent to pca.reduce().rows()
  size_type rank() const { return m_reduce.rows(); }

  /// \brief Orthogonal transformation matrix from a point in full [input...]
  ///        space to dimension-reduced [range(input...)] space.
  ///
  /// - Example use:
  /// \code
  /// Eigen::VectorXd reduced_point = pca.reduce()*full_point;
  /// \endcode
  ///
  Eigen::MatrixXd reduce() const { return m_reduce; }

  /// \brief Orthogonal transformation matrix from a point in dimension-reduced
  ///        [range(input...)] space to full [input...] space.
  ///
  /// - Example use:
  /// \code
  /// Eigen::VectorXd full_point = pca.expand()*reduced_point;
  /// \endcode
  ///
  Eigen::MatrixXd expand() const { return m_reduce.transpose(); }

 private:
  // transform full dimension [input...] column vector onto subspace
  // [range(input...)]
  Eigen::MatrixXd m_reduce;
};

/// \brief Construct a matrix consisting of blocks M and Identity(n,n)
///
/// \code
/// M ->  [ M 0 ]
///       [ 0 I ]
/// \endcode
inline Eigen::MatrixXd pad(const Eigen::MatrixXd &M, int n) {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(M.rows() + n, M.cols() + n);
  result << M, Eigen::MatrixXd::Zero(M.rows(), n),
      Eigen::MatrixXd::Zero(n, M.cols()), Eigen::MatrixXd::Identity(n, n);
  return result;
}

}  // namespace CASM

#endif
