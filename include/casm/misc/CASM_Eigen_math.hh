#ifndef CASM_Eigen_math
#define CASM_Eigen_math

#include <vector>
#include <cmath>
#include <complex>
#include <cassert>
#include <type_traits>
#include <boost/math/special_functions/round.hpp>

#include "./KroneckerTensorProduct.h"
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"
#include "casm/misc/CASM_math.hh"

#include "casm/external/Eigen/Eigenvalues"

namespace Local {
  template<typename T>
  inline static bool _compare(T const &a, T const &b, double tol) {
    return (a + tol) < b;
  }

  template<>
  inline bool _compare(std::complex<double> const &a, std::complex<double> const &b, double tol) {
    if(!CASM::almost_equal(real(a), real(b), tol)) {
      return real(a) < real(b);
    }
    else if(!CASM::almost_equal(imag(a), imag(b), tol)) {
      return imag(a) < imag(b);
    }
    return false;
  }
}
namespace CASM {

  /// \brief Floating point lexicographical comparison with tol
  inline bool float_lexicographical_compare(const Eigen::Ref<const Eigen::MatrixXd> &A, const Eigen::Ref<const Eigen::MatrixXd> &B, double tol) {
    return float_lexicographical_compare(A.data(), A.data() + A.size(), B.data(), B.data() + B.size(), tol);
  }

  /// \brief Floating point lexicographical comparison with tol
  template<typename Derived>
  bool colmajor_lex_compare(const Eigen::MatrixBase<Derived> &A, const Eigen::MatrixBase<Derived> &B, double tol) {
    if(A.cols() == B.cols()) {
      if(A.rows() == B.rows()) {
        for(Index j = 0; j < A.cols(); ++j) {
          for(Index i = 0; i < A.rows(); ++i) {
            if(!almost_equal(A(i, j), B(i, j), tol))
              return Local::_compare(A(i, j), B(i, j), tol);
          }
        }
        return false;
      }
      else if(A.rows() < B.rows())
        return true;
      else
        return false;
    }
    else if(A.cols() < B.cols())
      return true;
    else
      return false;
  }


  /// Returns reduced_column_echelon form of M. @param _tol is used to identify zero values.
  template<typename Derived>
  typename Derived::PlainObject reduced_column_echelon(Eigen::MatrixBase<Derived> const &M, double _tol) {
    typename Derived::PlainObject R(M);
    Index col = 0;
    for(Index row = 0; row < R.rows(); ++row) {
      Index i = 0;
      for(i = col; i < R.cols(); ++i) {
        if(!almost_zero(R(row, i), _tol)) {
          if(i != col)
            R.col(i).swap(R.col(col));
          break;
        }
      }
      if(i == R.cols())
        continue;
      R.col(col) /= R(row, col);
      for(i = 0; i < R.cols(); ++i) {
        if(i != col)
          R.col(i) -= R(row, i) * R.col(col);
      }
      ++col;
    }
    return R.leftCols(col);
  }


  Eigen::VectorXd eigen_vector_from_string(const std::string &tstr, const int &size);

  // *******************************************************************************************

  // Finds optimal assignments, based on cost_matrix, and returns total optimal cost
  double hungarian_method(const Eigen::MatrixXd &cost_matrix, std::vector<Index> &optimal_assignments, const double _tol);


  ///\brief Return pointer one past end of vector. Equivalent to convainer.data()+container.size()
  template<typename Derived>
  typename Derived::Scalar *end_ptr(Eigen::PlainObjectBase<Derived> &container) {
    return container.data() + container.size();
  }

  ///\brief Return const pointer one past end of vector. Equivalent to convainer.data()+container.size()
  template<typename Derived>
  typename Derived::Scalar const *end_ptr(Eigen::PlainObjectBase<Derived> const &container) {
    return container.data() + container.size();
  }

  /** \defgroup LinearAlgebra
  *
  *  \ingroup Container
  *  \brief Linear algebra routines
  *
  *  @{
  */

  void get_Hermitian(Eigen::MatrixXcd &original_mat, Eigen::MatrixXcd &hermitian_mat, Eigen::MatrixXcd &antihermitian_mat); //Ivy
  bool is_Hermitian(Eigen::MatrixXcd &mat); //Ivy

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
  typename Derived::Scalar triple_product(const Derived &vec0,
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

  /// \brief Round Eigen::MatrixXd
  ///
  /// \returns an Eigen:MatrixXd
  ///
  /// \param M Eigen::MatrixXd to be rounded
  ///
  /// For each coefficient, sets \code M(i,j) = boost::math::round(Mdouble(i, j)) \endcode
  ///
  template<typename Derived>
  Eigen::CwiseUnaryOp< decltype(std::ptr_fun(boost::math::round<typename Derived::Scalar>)), const Derived >
  round(const Eigen::MatrixBase<Derived> &val) {
    return val.unaryExpr(std::ptr_fun(boost::math::round<typename Derived::Scalar>));
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
  Eigen::CwiseUnaryOp< decltype(std::ptr_fun(boost::math::iround<typename Derived::Scalar>)), const Derived >
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
  Eigen::CwiseUnaryOp< decltype(std::ptr_fun(boost::math::lround<typename Derived::Scalar>)), const Derived >
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
    int i, j;
    for(i = 0; i < row; ++i) {
      for(j = 0; j < col; ++j) {
        subM(i, j) = M(i, j);
      }
      for(++j; j < M.cols(); ++j) {
        subM(i, j - 1) = M(i, j);
      }
    }

    for(++i; i < M.rows(); ++i) {
      for(j = 0; j < col; ++j) {
        subM(i - 1, j) = M(i, j);
      }
      for(++j; j < M.cols(); ++j) {
        subM(i - 1, j - 1) = M(i, j);
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
        Mcof(i, j) = ((((i + j) % 2) == 0) ? matrix_minor(M, i, j) : -matrix_minor(M, i, j));
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

  Eigen::Matrix3d polar_decomposition(Eigen::Matrix3d const &mat);

  std::vector<Eigen::Matrix3i> _unimodular_matrices(bool positive, bool negative, int range = 1);
  const std::vector<Eigen::Matrix3i> &positive_unimodular_matrices();
  const std::vector<Eigen::Matrix3i> &negative_unimodular_matrices();

  template<int range = 1>
  const std::vector<Eigen::Matrix3i> &unimodular_matrices() {
    static std::vector<Eigen::Matrix3i> static_all(_unimodular_matrices(true, true, range));
    return static_all;
  }


  /** @} */

}

namespace Eigen {
  /// \brief Round Eigen::MatrixXd
  ///
  /// \returns an Eigen:MatrixXd
  ///
  /// \param M Eigen::MatrixXd to be rounded
  ///
  /// For each coefficient, sets \code M(i,j) = boost::math::floor(Mdouble(i, j)) \endcode
  ///
  /* namespace Local { */
  /*   template<typename T> */
  /*   struct _Floor { */
  /*     T operator()(T val)const { */
  /*       return floor(val); */
  /*     } */
  /*   }; */
  /* } */
  /* template<typename Derived> */
  /* CwiseUnaryOp<Local::_Floor<typename Derived::Scalar>, const Derived > */
  /* floor(const MatrixBase<Derived> &val) { */
  /*   return val.unaryExpr(Local::_Floor<typename Derived::Scalar>()); */
  /* } */


  /// \brief Equivalent to almost_zero(double(val.norm()), tol);
  template <typename Derived>
  bool almost_zero(const MatrixBase<Derived> &val, double tol = CASM::TOL) {
    return val.isZero(tol);
  }

  template<typename Derived>
  double length(const Eigen::MatrixBase<Derived> &value) {
    return value.norm();
  }

  //*******************************************************************************************
  ///Take a vector of doubles, and multiply by some factor that turns it into a vector of integers (within a tolerance)
  template <typename Derived>
  Eigen::Matrix<int,
        Derived::RowsAtCompileTime,
        Derived::ColsAtCompileTime>
  scale_to_int(const Eigen::MatrixBase<Derived> &val, double _tol = CASM::TOL) {

    using CASM::Index;
    using CASM::almost_zero;

    typedef Eigen::Matrix<int, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime> int_mat_type;
    typedef Eigen::Matrix<double, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime> dub_mat_type;

    int_mat_type ints(int_mat_type::Zero(val.rows(), val.cols()));

    dub_mat_type dubs(val);

    Index min_i(-1);
    double min_coeff = 2; //all values are <=1;
    for(Index i = 0; i < dubs.rows(); i++) {
      for(Index j = 0; j < dubs.cols(); j++) {
        if(almost_zero(dubs(i, j))) {
          dubs(i, j) = 0.0;
        }
        else if(std::abs(dubs(i, j)) < std::abs(min_coeff)) {
          min_coeff = dubs(i, j);
          min_i = i;
        }
      }
    }
    if(CASM::valid_index(min_i))
      dubs /= std::abs(min_coeff);
    else
      return ints;


    //We want to multiply the miller indeces by some factor such that all indeces become integers.
    //In order to do this we pick a tolerance to work with and round the miller indeces if they are close
    //enough to the integer value (e.g. 2.95 becomes 3). Choosing a tolerance that is too small will
    //result in the "primitive-slab" blowing up.

    //Begin choosing a factor and multiply all indeces by it (starting with 1). Then round the non-smallest
    //miller indeces (smallest index requires no rounding, since it will always be a perfect
    //integer thanks to the previous division).
    //Next take absolute value of difference between rounded indeces and actual values (int_diff 1 & 2).
    //If the difference for both indeces is smaller than the tolerance then you've reached the desired
    //accuracy and the rounded indeces can be used to construct the "primitive-slab" cell. If not, increase the
    //factor by 1 and try again, until the tolerance is met.
    bool within_tol = false;

    dub_mat_type tdubs;
    Index i, j;
    for(Index factor = 1; factor < 1000 && !within_tol; factor++) {
      tdubs = double(factor) * dubs;
      for(i = 0; i < dubs.rows(); i++) {
        for(j = 0; j < dubs.cols(); j++) {

          if(!CASM::almost_zero(boost::math::round(tdubs(i, j)) - tdubs(i, j), _tol))
            break;
        }
        if(j < dubs.cols())
          break;
      }
      if(dubs.rows() <= i)
        within_tol = true;
    }

    if(within_tol) {
      for(i = 0; i < dubs.rows(); i++) {
        for(j = 0; j < dubs.cols(); j++) {
          ints(i, j) = boost::math::round(tdubs(i, j));
        }
      }
    }

    return ints;
  }

  template <typename Derived1, typename Derived2>
  inline
  bool almost_equal(const Eigen::MatrixBase<Derived1> &val1, const Eigen::MatrixBase<Derived2> &val2, double tol = CASM::TOL) {
    return almost_zero(val1 - val2, tol);
  }

  /**
   * Checks to see whether the given matrix is symmetric
   * by checking if its transpose is equal to itself.
   * Only works for square matrices n x n.
   * (Reflected along 0,0 to n,n)
   */

  template <typename Derived>
  inline
  bool is_symmetric(const Eigen::MatrixBase<Derived> &test_mat, double test_tol = CASM::TOL) {
    return almost_zero(test_mat - test_mat.transpose(), test_tol);
  }

  /**
   * Checks to see if the given matrix is persymmetric, i.e.
   * whether it's symmetric along the cross diagonal.
   * Only works for square matrices n x n.
   * (Reflected along 0,n to n,0)
   */

  template <typename Derived>
  inline
  bool is_persymmetric(const Eigen::MatrixBase<Derived> &test_mat, double test_tol = CASM::TOL) {
    //Reverse order of columns and rows
    auto rev_mat = test_mat.colwise().reverse().eval().rowwise().reverse().eval();
    return almost_zero(test_mat - rev_mat.transpose(), test_tol);
  }

  /**
   * Checks to see if the given matrix is bisymmetric, i.e.
   * whether it's symmetric along both diagonals.
   * Only works for square matrices n x n.
   * (Reflected along 0,n to n,0 AND 0,0 to n,n)
   */

  template <typename Derived>
  inline
  bool is_bisymmetric(const Eigen::MatrixBase<Derived> &test_mat, double test_tol = CASM::TOL) {
    return (is_symmetric(test_mat, test_tol) && is_persymmetric(test_mat, test_tol));
  }

  template<typename Derived>
  std::vector<CASM::Index> partition_distinct_values(const Eigen::MatrixBase<Derived> &value, double tol = CASM::TOL) {
    std::vector<CASM::Index> subspace_dims;
    CASM::Index last_i = 0;
    for(CASM::Index i = 1; i < value.size(); i++) {
      if(!CASM::almost_equal(value[last_i], value[i], tol)) {
        subspace_dims.push_back(i - last_i);
        last_i = i;
      }
    }
    subspace_dims.push_back(value.size() - last_i);
    return subspace_dims;

  }

  template<typename Derived, typename Index = typename Derived::Index>
  Index print_matrix_width(std::ostream &s, const Derived &m, Index width) {
    for(Index j = 0; j < m.cols(); ++j) {
      for(Index i = 0; i < m.rows(); ++i) {
        std::stringstream sstr;
        sstr.copyfmt(s);
        sstr << m.coeff(i, j);
        width = std::max<Index>(width, Index(sstr.str().length()));
      }
    }
    return width;
  }


}

#endif
