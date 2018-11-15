#ifndef CASM_Eigen_math
#define CASM_Eigen_math

#include <vector>
#include "casm/CASM_global_Eigen.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {

  /// \brief Floating point lexicographical comparison with tol
  inline bool float_lexicographical_compare(const Eigen::VectorXd &A, const Eigen::VectorXd &B, double tol) {
    return float_lexicographical_compare(A.data(), A.data() + A.size(), B.data(), B.data() + B.size(), tol);
  }


  Eigen::VectorXd eigen_vector_from_string(const std::string &tstr, const int &size);

  // *******************************************************************************************

  void poly_fit(Eigen::VectorXcd &xvec, Eigen::VectorXcd &yvec, Eigen::VectorXcd &coeffs, int degree); //Ivy

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

}

namespace Eigen {

  /// \brief Equivalent to almost_zero(double(val.norm()), tol);
  template <typename Derived>
  bool almost_zero(const Eigen::MatrixBase<Derived> &val, double tol = CASM::TOL) {
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

    typedef Eigen::Matrix<int, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime> int_mat_type;
    typedef Eigen::Matrix<double, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime> dub_mat_type;

    int_mat_type ints(int_mat_type::Zero(val.rows(), val.cols()));

    dub_mat_type dubs(val);

    Index min_i(-1), min_j(-1);
    double min_coeff = 2; //all values are <=1;
    for(Index i = 0; i < dubs.rows(); i++) {
      for(Index j = 0; j < dubs.cols(); j++) {
        if(CASM::almost_zero(dubs(i, j))) {
          dubs(i, j) = 0.0;
        }
        else if(std::abs(dubs(i, j)) < std::abs(min_coeff)) {
          min_coeff = dubs(i, j);
          min_i = i;
          min_j = j;
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
          if(!CASM::almost_zero(round(tdubs(i, j)) - tdubs(i, j), _tol))
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
          ints(i, j) = round(tdubs(i, j));
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
}

#endif
