//
//  LinearAlgebra.hh
//  CASM
//
//

#ifndef LINEARALGEBRA_HH
#define LINEARALGEBRA_HH

#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <cassert>

#include <boost/math/special_functions/round.hpp>

#include "casm/external/Eigen/Dense"
#include "casm/external/Eigen/Eigenvalues"

#include "casm/container/Array.hh"
#include "casm/misc/CASM_math.hh"
namespace CASM {

  /** \defgroup LinearAlgebra
   *
   *  \ingroup Container
   *  \brief Linear algebra routines
   *
   *  @{
  */

  void get_Hermitian(Eigen::MatrixXcd &original_mat, Eigen::MatrixXcd &hermitian_mat, Eigen::MatrixXcd &antihermitian_mat); //Ivy
  bool is_Hermitian(Eigen::MatrixXcd &mat); //Ivy
  void poly_fit(Eigen::VectorXcd &xvec, Eigen::VectorXcd &yvec, Eigen::VectorXcd &coeffs, int degree);

  /// \brief Return the hermite normal form, M == H*V
  std::pair<Eigen::MatrixXi, Eigen::MatrixXi> hermite_normal_form(const Eigen::MatrixXi &M);

  /// \brief Get angle, in radians, between two vectors on range [0,pi]
  double angle(const Eigen::Ref<const Eigen::Vector3d> &a, const Eigen::Ref<const Eigen::Vector3d> &b);

  /// \brief signed angle, in radians, between -pi and pi that describe separation in direction of two vectors
  double signed_angle(const Eigen::Ref<const Eigen::Vector3d> &a,
                      const Eigen::Ref<const Eigen::Vector3d> &b,
                      const Eigen::Ref<const Eigen::Vector3d> &pos_ref);

  /// \brief Round entries that are within tol of being integer to that integer value
  Eigen::MatrixXd pretty(const Eigen::MatrixXd &M, double tol);


  template<typename Derived>
  typename Derived::Scalar triple_prod(const Derived &vec0,
                                       const Derived &vec1,
                                       const Derived &vec2) {
    return (vec0.cross(vec1)).dot(vec2);
  }

  /// \brief Check if Eigen::Matrix is integer
  template<typename Derived>
  bool is_integer(const Eigen::MatrixBase<Derived> &M, double tol) {
    for(Index i = 0; i < M.rows(); i++) {
      for(Index j = 0; j < M.cols(); j++) {
        if(!almost_zero(boost::math::lround(M(i, j)) - M(i, j), tol))
          return false;
      }
    }
    return true;
  }

  /// \brief Check if Eigen::Matrix is unimodular
  template<typename Derived>
  bool is_unimodular(const Eigen::MatrixBase<Derived> &M, double tol) {
    return is_integer(M, tol) && almost_equal(std::abs(M.determinant()), static_cast<typename Derived::Scalar>(1), tol);
  }

  /// \brief Check if Eigen::Matrix is diagonal
  template<typename Derived>
  bool is_diagonal(const Eigen::MatrixBase<Derived> &M, double tol = TOL) {
    return M.isDiagonal(tol);
  }

  /// \brief Round Eigen::MatrixXd to Eigen::MatrixXi
  ///
  /// \returns an Eigen:MatrixXi
  ///
  /// \param M Eigen::MatrixXd to be rounded to integer
  ///
  /// For each coefficient, sets \code Mint(i,j) = boost::math::iround(Mdouble(i, j)) \endcode
  ///
  template<typename Derived>
  Eigen::CwiseUnaryOp< decltype(std::ptr_fun(boost::math::iround<typename Derived::Scalar>)) , const Derived >
  iround(const Eigen::MatrixBase<Derived> &val) {
    return val.unaryExpr(std::ptr_fun(boost::math::iround<typename Derived::Scalar>));
  }

  /// \brief Round Eigen::MatrixXd to Eigen::MatrixXl
  ///
  /// \returns an Eigen:MatrixXl
  ///
  /// \param M Eigen::MatrixXd to be rounded to integer
  ///
  /// For each coefficient, sets \code Mint(i,j) = std::lround(Mdouble(i, j)) \endcode
  ///
  template<typename Derived>
  Eigen::CwiseUnaryOp< decltype(std::ptr_fun(boost::math::lround<typename Derived::Scalar>)) , const Derived >
  lround(const Eigen::MatrixBase<Derived> &val) {
    return val.unaryExpr(std::ptr_fun(boost::math::lround<typename Derived::Scalar>));
  }

  /// \brief Return the minor of integer Matrix M element row, col
  ///
  /// \returns the minor of element row, col
  ///
  /// \param M matrix
  /// \param row,col element row and column
  ///
  /// The minor of element row, col of a matrix is the determinant of the submatrix of
  /// M which does not include any elements from row 'row' or column 'col'
  ///
  template<typename Derived>
  typename Derived::Scalar matrix_minor(const Eigen::MatrixBase<Derived> &M, int row, int col) {

    // create the submatrix of M which doesn't include any elements from M in 'row' or 'col'
    Eigen::Matrix < typename Derived::Scalar,
          Derived::RowsAtCompileTime - 1,
          Derived::ColsAtCompileTime - 1 >
          subM(M.rows() - 1, M.cols() - 1);
    int _i, _j;
    for(int i = 0; i < M.rows(); i++) {
      for(int j = 0; j < M.cols(); j++) {
        if(i == row || j == col)
          continue;
        (i < row) ? _i = i : _i = i - 1;
        (j < col) ? _j = j : _j = j - 1;
        subM(_i, _j) = M(i, j);
      }
    }

    // return determinant of subM
    return subM.determinant();
  }

  /// \brief Return cofactor matrix
  ///
  /// \returns cofactor matrix
  ///
  /// \param M matrix
  ///
  /// The cofactor matrix is \code C(i,j) = pow(-1,i+j)*matrix_minor(M,i,j) \endcode
  ///
  template<typename Derived>
  Eigen::Matrix<typename Derived::Scalar,
        Derived::RowsAtCompileTime,
        Derived::ColsAtCompileTime>
  cofactor(const Eigen::MatrixBase<Derived> &M) {
    Eigen::Matrix<typename Derived::Scalar,
          Derived::RowsAtCompileTime,
          Derived::ColsAtCompileTime>
          Mcof(M.rows(), M.cols());
    for(int i = 0; i < Mcof.rows(); i++) {
      for(int j = 0; j < Mcof.cols(); j++) {
        (((i + j) % 2) == 0) ? Mcof(i, j) = matrix_minor(M, i, j) : Mcof(i, j) = -matrix_minor(M, i, j);
      }
    }
    return Mcof;
  }

  /// \brief Return adjugate matrix
  ///
  /// \returns adjugate matrix
  ///
  /// \param M matrix
  ///
  /// The adjugate matrix is the transpose of the cofactor matrix.
  /// The adjugate matrix is related to the inverse of a matrix through
  /// \code M.inverse() == adjugate(M)/M.determinant() \endcode
  ///
  // Note: If Eigen ever implements integer factorizations, we should probably change how this works
  template<typename Derived>
  Eigen::Matrix<typename Derived::Scalar,
        Derived::RowsAtCompileTime,
        Derived::ColsAtCompileTime>
  adjugate(const Eigen::MatrixBase<Derived> &M) {
    return cofactor(M).transpose();
  }

  namespace normal_form_impl {
    /// \brief Return an elementary hermite transformation
    ///
    /// \returns an integer matrix, E, with determinant 1 such that E * [a;b] = [g;0]
    ///
    /// \param a,b integers
    ///
    /// \see smith_normal_form
    ///
    template<typename Scalar>
    Eigen::Matrix<Scalar, 3, 3> _elementary_hermite_op(Scalar a, Scalar b, Scalar i, Scalar j) {
      typedef Eigen::Matrix<Scalar, 3, 3> matrix_type;
      matrix_type tmat = matrix_type::Identity();
      Scalar tgcf = extended_gcf(a, b, tmat(i, i), tmat(i, j));
      if(tgcf == 0) {
        return matrix_type::Identity();
      }
      tmat(j, i) = -b / tgcf;
      tmat(j, j) = a / tgcf;
      return tmat;
    }

  }

  /// \brief Return the integer inverse matrix of an invertible integer matrix
  ///
  /// \returns a 3x3 integer matrix inverse
  ///
  /// \param M an 3x3 invertible integer matrix
  ///
  template<typename Derived>
  Eigen::Matrix<typename Derived::Scalar,
        Derived::RowsAtCompileTime,
        Derived::ColsAtCompileTime>
  inverse(const Eigen::MatrixBase<Derived> &M) {
    return adjugate(M) / M.determinant();
  }

  /// \brief Return the smith normal form, M == U*S*V
  ///
  /// \param M 3x3 integer matrix
  /// \param[out] U set to U: is integer, U.determinant() == 1
  /// \param[out] S set to S: is diagonal, integer, and nonnegative, S(i,i) divides S(i+1,i+1) for all i
  /// \param[out] V set to V: is integer, V.determinant() == 1
  ///
  /// Adapted from Matlab implementation written by John Gilbert (gilbert@parc.xerox.com):
  /// - http://www.mathworks.com/matlabcentral/newsreader/view_thread/13728
  ///
  template<typename DerivedIn, typename DerivedOut>
  void smith_normal_form(const Eigen::MatrixBase<DerivedIn> &M,
                         Eigen::MatrixBase<DerivedOut> &U,
                         Eigen::MatrixBase<DerivedOut> &S,
                         Eigen::MatrixBase<DerivedOut> &V) {

    static_assert(std::is_same<typename DerivedIn::Scalar, typename DerivedOut::Scalar>::value,
                  "ALL ARGUMENTS TO CASM::smith_normal_form() MUST BE BASED ON SAME SCALAR TYPE!");
    using namespace normal_form_impl;
    typedef typename DerivedOut::Scalar Scalar;

    U = V = DerivedOut::Identity();
    S = M;
    //std::cout << "S is:\n" << S << "\n\n";

    DerivedOut tmat = U;

    int i, j;

    // Bidiagonalize S with elementary Hermite transforms.
    for(j = 0; j < 3; j++) {
      //take j as column index and zero elements of j below the diagonal
      for(i = j + 1; i < 3; i++) {
        if(S(i, j) == 0) continue;
        tmat = _elementary_hermite_op<Scalar>(S(j, j), S(i, j), j, i);
        S = tmat * S;
        U = U * inverse(tmat);
      }

      //take j as row index and zero elements of j after superdiagonal
      for(i = j + 2; i < 3; i++) {
        if(S(j, i) == 0) continue;
        tmat = _elementary_hermite_op<Scalar>(S(j, j + 1), S(j, i), j + 1, i);
        S = S * tmat.transpose();
        V = inverse(tmat.transpose()) * V;
      }
    }

    // std::cout << "About to hunt zeros, SNF decomposition is currently:\n" << U << "\n\n" << S << "\n\n" << V << "\n\n";
    while(!S.isDiagonal()) {
      //find off-diagonal element
      int b = 0;
      for(b = 0; b < 2; b++) {
        if(S(b, b + 1))
          break;
      }

      // To guarantee reduction in S(b,b), first make S(b,b) positive
      // and make S(b,b+1) nonnegative and less than S(b,b).
      if(S(b, b) < 0) {
        for(i = 0; i < 3; i++) {
          S(b, i) = -S(b, i);
          U(i, b) = -U(i, b);
        }
      }

      //std::cout << "S before q:\n" << S << '\n';

      if(S(b, b)) {
        Scalar q = S(b, b + 1) / S(b, b);
        if(S(b, b + 1) % S(b, b) < 0)
          q -= 1;
        tmat = DerivedOut::Identity();
        tmat(b + 1, b) = -q;
        S = S * tmat.transpose();
        V = inverse(tmat.transpose()) * V;
        //std::cout << "tmat for q:\n" << tmat << '\n';
        //std::cout << "S after q:\n" << S << '\n';
      }
      else {
        tmat = DerivedOut::Identity();
        tmat(b, b) = 0;
        tmat(b, b + 1) = 1;
        tmat(b + 1, b + 1) = 0;
        tmat(b + 1, b) = 1;
        S = S * tmat.transpose();
        V = inverse(tmat.transpose()) * V;

      }

      if(!S(b, b + 1))
        continue;

      tmat = _elementary_hermite_op<Scalar>(S(b, b), S(b, b + 1), b, b + 1);
      S = S * tmat.transpose();
      V = inverse(tmat.transpose()) * V;
      for(j = 0; j < 2; j++) {
        tmat = _elementary_hermite_op<Scalar>(S(j, j), S(j + 1, j), j, j + 1);
        S = tmat * S;
        U = U * inverse(tmat);

        if(j + 2 >= 3)
          continue;

        tmat = _elementary_hermite_op<Scalar>(S(j, j + 1), S(j, j + 2), j + 1, j + 2);
        S = S * tmat.transpose();
        V = inverse(tmat.transpose()) * V;
      }
    }

    //Now it's diagonal -- make it non-negative
    for(j = 0; j < 3; j++) {
      if(S(j, j) < 0) {
        for(i = 0; i < 3; i++) {
          S(j, i) = -S(j, i);
          U(i, j) = -U(i, j);
        }
      }
    }

    //sort diagonal elements
    for(i = 0; i < 3; i++) {
      for(j = i + 1; j < 3; j++) {
        if(S(i, i) > S(j, j)) {
          S.row(i).swap(S.row(j));
          S.col(i).swap(S.col(j));
          U.col(i).swap(U.col(j));
          V.row(i).swap(V.row(j));
        }
      }
    }
    //enforce divisibility condition:
    for(i = 0; i < 3; i++) {
      if(S(i, i) == 0) continue;
      for(j = i + 1; j < 3; j++) {
        if(S(j, j) % S(i, i) == 0) continue;
        //Replace S(i,i), S(j,j) by their gcd and lcm respectively.
        tmat = DerivedOut::Identity();
        DerivedOut tmat2(tmat);
        Scalar a(S(i, i)), b(S(j, j)), c, d, tgcf;
        tgcf = extended_gcf(a, b, c, d);
        tmat(i, i) = 1;
        tmat(i, j) = d;
        tmat(j, i) = -b / tgcf;
        tmat(j, j) = (a * c) / tgcf;

        tmat2(i, i) = c;
        tmat2(i, j) = 1;
        tmat2(j, i) = -(b * d) / tgcf;
        tmat2(j, j) = a / tgcf;
        S = tmat * S * tmat2.transpose();
        U = U * inverse(tmat);
        V = inverse(tmat2.transpose()) * V;
      }
    }
    return;
  }

  std::vector<Eigen::Matrix3i> _unimodular_matrices(bool positive, bool negative);
  const std::vector<Eigen::Matrix3i> &positive_unimodular_matrices();
  const std::vector<Eigen::Matrix3i> &negative_unimodular_matrices();
  const std::vector<Eigen::Matrix3i> &unimodular_matrices();

  /** @} */
}
#endif
