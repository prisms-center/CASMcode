#ifndef CASM_MATH_HH
#define CASM_MATH_HH
#include "casm/CASM_global_definitions.hh"
#include "casm/container/Array.hh"
//Maybe we should transition to boost math library?
//#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <iostream>
#include <cmath>
#include <cstddef>
#include <complex>
#include <string>
#include <sstream>

namespace CASM {

  // *******************************************************************************************

  //round function -- rounds to nearest whole number; half-integers are rounded away from zero
  int round(double val);

  // *******************************************************************************************

  template <typename T>
  T min(const T &A, const T &B) {
    return A < B ? A : B;
  }

  // *******************************************************************************************

  template <typename T>
  T max(const T &A, const T &B) {
    return A > B ? A : B;
  }

  // *******************************************************************************************

  /// \brief If T is not integral, use std::abs(val) < tol;
  template <typename T, typename std::enable_if<std::is_floating_point<T>::value, T>::type * = nullptr>
  inline
  bool almost_zero(const T &val, double tol = TOL) {
    return std::abs(val) < tol;
  }

  /// \brief If std::complex<T>, use std::abs(val) < tol;
  template <typename T, typename std::enable_if<std::is_floating_point<T>::value, T>::type * = nullptr>
  inline
  bool almost_zero(const std::complex<T> &val, double tol = TOL) {
    return std::abs(val) < tol;
  }

  /// \brief If T is integral, val == 0;
  template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
  inline
  bool almost_zero(const T &val, double tol = TOL) {
    return val == 0;
  }

  /// \brief Equivalent to almost_zero(double(val.norm()), tol);
  template <typename Derived>
  bool almost_zero(const Eigen::MatrixBase<Derived> &val, double tol = TOL) {
    return val.isZero(tol);
  }

  // *******************************************************************************************

  /// \brief If T is not integral, use almost_zero(val1 - val2, tol);
  template < typename T, typename std::enable_if < !std::is_integral<T>::value, T >::type * = nullptr >
  bool almost_equal(const T &val1, const T &val2, double tol = TOL) {
    return almost_zero(val1 - val2, tol);
  }

  /// \brief If T is integral type, use val1 == val2;
  template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
  bool almost_equal(const T &val1, const T &val2, double tol = TOL) {
    return val1 == val2;
  }

  // *******************************************************************************************

  /// \brief Floating point comparison with tol, return A < B
  ///
  /// Implements:
  /// \code
  /// if(!almost_equal(A,B,tol)) { return A < B; }
  /// return false;
  /// \endcode
  template <typename T, typename std::enable_if<std::is_floating_point<T>::value, T>::type * = nullptr>
  bool compare(const T &A, const T &B, double tol) {
    if(!almost_equal(A, B, tol)) {
      return A < B;
    }
    return false;
  }

  struct FloatCompare {

    double tol;

    FloatCompare(double _tol) : tol(_tol) {}

    template<typename T>
    bool operator()(const T &A, const T &B) const {
      return compare(A, B, tol);
    }
  };


  /// \brief Floating point lexicographical comparison with tol
  template<class InputIt1, class InputIt2>
  bool float_lexicographical_compare(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, double tol) {
    FloatCompare compare(tol);
    return std::lexicographical_compare(first1, last1, first2, last2, compare);
  }

  /// \brief Floating point lexicographical comparison with tol
  inline bool float_lexicographical_compare(const Eigen::VectorXd &A, const Eigen::VectorXd &B, double tol) {
    return float_lexicographical_compare(A.data(), A.data() + A.size(), B.data(), B.data() + B.size(), tol);
  }

  // *******************************************************************************************
  //Return sign of number

  template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
  int sgn(T val) {
    return (T(0) < val) - (val < T(0));
  }

  template <typename T, typename std::enable_if<std::is_floating_point<T>::value, T>::type * = nullptr>
  int float_sgn(T val, double compare_tol = TOL) {
    T zeroval(0);
    if(compare(zeroval, val, compare_tol)) {
      return 1;
    }

    else if(compare(val, zeroval, compare_tol)) {
      return -1;
    }

    else {
      return 0;
    }
  }
  // *******************************************************************************************

  // Works for signed and unsigned types
  template <typename IntType>
  IntType nchoosek(IntType n, IntType k) {
    assert(k <= n && 0 <= k);
    if(n < 2 * k)
      k = n - k;

    n -= k;

    IntType result(1);
    for(IntType i = 1; i < k + 1; i++) {
      result *= n + i;
      result /= i;
    }

    return result;
  }


  // *******************************************************************************************
  /// \brief Computes the Damerescau-Levenshtein distance
  ///   -- the number of edits (deletions, insertions, transpositions) to go from string 'a' to string 'b'
  int dl_string_dist(const std::string &a, const std::string &b);

  // *******************************************************************************************

  double ran0(int &idum);

  // *******************************************************************************************

  using boost::math::factorial;

  // *******************************************************************************************
  /// Find greatest common factor
  int gcf(int i1, int i2);

  // *******************************************************************************************
  /// Find least common multiple
  int lcm(int i1, int i2);

  // *******************************************************************************************

  // evaluates Gaussians using the formula:
  // f(x) = a*e^(-(x-b)^2/(c^2))
  double gaussian(double a, double x, double b, double c);

  // calculates Gaussian moments given by the integral:
  // m = \int_{\infty}^{\infty} dx x^pow*exp[-x^2/(2*sigma^2)]/(\sqrt(2*\pi)*sigma)
  double gaussian_moment(int expon, double sigma);

  // calculates Gaussian moments given by the integral:
  // m = \int_{\infty}^{\infty} dx x^pow*exp[-(x-x0)^2/(2*sigma^2)]/(\sqrt(2*\pi)*sigma)
  double gaussian_moment(int expon, double sigma, double x0);


  Eigen::VectorXd eigen_vector_from_string(const std::string &tstr, const int &size);

  // *******************************************************************************************
  /// \brief Calculate greatest common factor of two integers, and bezout coefficients
  ///
  /// \returns greatest common factor
  ////
  /// \param i1,i2 two integers for which to find greatest common factor
  /// \param[out] p1, p2 bezout coefficients such that p1*i1 + p2*i2 = gcf(abs(i1),abs(i2));
  ///
  /// \see smith_normal_form
  ///
  template<typename IntType>
  IntType extended_gcf(IntType i1, IntType i2, IntType &p1, IntType &p2) {
    IntType s1 = sgn(i1);
    IntType s2 = sgn(i2);

    IntType a = i1, b = i2;
    p1 = 1;
    p2 = 0;
    IntType tp1(0), tp2(1), quotient;

    i1 = std::abs(i1);
    i2 = std::abs(i2);

    while(i2 != 0) {

      quotient = i1 / i2;
      i1 = i1 % i2;
      swap(i1, i2);
      p1 = p1 - quotient * tp1;
      swap(p1, tp1);
      p2 = p2 - quotient * tp2;
      swap(p2, tp2);
    }

    p1 *= s1;
    p2 *= s2;
    return p1 * a + p2 * b;
  }

  // *******************************************************************************************
  /// Find least common multiple
  int lcm(int i1, int i2);


  // *******************************************************************************************
  // finds rational number that approximates 'val' to within 'tol' --
  //      --  almost_int(double(denominator)*val, tol) == true OR almost_int(val/double(numerator), tol) == true OR
  //      --  denominator is always positive
  //      --  sign(numerator)=sign(val)
  //      --  if(almost_zero(val,tol)) --> numerator = 0, denominator = 1
  // *******************************************************************************************

  void nearest_rational_number(double val, long &numerator, long &denominator, double tol = TOL);

  // *******************************************************************************************
  /*
   * Finds best irrational number approximation of double 'val' and
   * returns tex-formated string that contains irrational approximation
   * Searches numbers of the form (x/y)^(1/z), where x and y range from 1 to 'lim'
   * z ranges from 1 to 'max_pow'
   */
  // *******************************************************************************************

  std::string irrational_to_tex_string(double val, int lim, int max_pow = 2);

  // *******************************************************************************************
  //John G 010413
  int mod(int a, int b);

  // *******************************************************************************************

  double cuberoot(double number);

  // *******************************************************************************************

  void poly_fit(Eigen::VectorXcd &xvec, Eigen::VectorXcd &yvec, Eigen::VectorXcd &coeffs, int degree); //Ivy

  // //////////////////////////////////////////
  // //////////////////////////////////////////
  // Array Function declarations:

  // *******************************************************************************************

  // returns 'i' if 'input' is equivalent to 'unique[i]', w.r.t. permutation of the equivalent elements of 'input'.
  // equivalent elements are specified by 'ind_equiv'
  // if 'input' specifies a new combination of integers, unique.size() is returned
  Index which_unique_combination(const Array<Index> &input, const Array<Index>::X2 &unique, const Array<Index>::X2 &ind_equiv);

  // *******************************************************************************************

  Index which_unique_combination(const Array<Index> &input, const Array<Index>::X2 &unique);

  // *******************************************************************************************
  /// Find least common multiple
  int lcm(const Array<int> &series);

  // *******************************************************************************************

  ReturnArray< Array<int> > get_prime_factors(int target);

  // *******************************************************************************************

  template <typename IntType>
  IntType multinomial_coeff(const Array<IntType> &exponents) {
    IntType tcoeff(1), tsum(0);
    for(Index i = 0; i < exponents.size(); i++) {
      tsum += exponents[i];
      tcoeff *= nchoosek(tsum, exponents[i]);
    }
    return tcoeff;
  }

  // *******************************************************************************************
  // get multinomial coefficient for only a subset of the coeffs
  template <typename IntType>
  IntType multinomial_coeff(const Array<IntType> &exponents, const Array<Index> &sublist) {
    IntType tcoeff(1), tsum(0);
    for(Index i = 0; i < sublist.size(); i++) {
      tsum += exponents[sublist[i]];
      tcoeff *= nchoosek(tsum, exponents[sublist[i]]);
    }
    return tcoeff;
  }
  // ************************************************************
  template <typename T>
  bool almost_equal(const Array<T> &A, const Array<T> &B, double tol = TOL) {
    if(A.size() != B.size())
      return false;
    for(Index i = 0; i < A.size(); i++) {
      if(!almost_equal(A[i], B[i], tol))
        return false;
    }
    return true;
  }

  // ************************************************************
  // cumulative sum, first element is 0 and final elment is array.sum()
  template <typename T>
  ReturnArray<T> cum_sum(const Array<T> &arr) {
    Array<T> result;
    result.reserve(arr.size() + 1);
    result.push_back(0);
    for(Index i = 0; i < arr.size(); i++) {
      result.push_back(result[i] + arr[i]);
    }
    return result;
  }

  // ************************************************************

  template<typename Derived>
  double length(const Eigen::MatrixBase<Derived> &value) {
    return value.norm();
  }

  // ************************************************************

  template<typename Derived>
  ReturnArray<Index> partition_distinct_values(const Eigen::MatrixBase<Derived> &value, double tol = TOL) {
    Array<Index> subspace_dims;
    Index last_i = 0;
    for(Index i = 1; i < value.size(); i++) {
      if(!almost_equal(value[last_i], value[i], tol)) {
        subspace_dims.push_back(i - last_i);
        last_i = i;
      }
    }
    subspace_dims.push_back(value.size() - last_i);
    return subspace_dims;
  }

  // Finds optimal assignments, based on cost_matrix, and returns total optimal cost
  double hungarian_method(const Eigen::MatrixXd &cost_matrix, std::vector<Index> &optimal_assignments, const double _tol);

  namespace HungarianMethod_impl {
    // *******************************************************************************************
    /* Hungarian Algorithm Routines
     * Step 1: reduce the rows by smallest element
     * Step 2: star zeros
     * Step 3: cover columns with starred zeros and check assignement
     *         if K columns covered DONE. Else goto 4.
     * Step 4: Find uncovered zero and prime. if no starred zero in row
     *         goto 5. Otherwise, cover row, uncover column with star zero.
     *         Continue until all zeros are covered. Store smalles
     *         uncovered value goto 6.
     * Step 5: Build alternating prime and star zeros. Goto 3.
     * Step 6: Add value from 4 to all covered rows, and subtract it
     *         from uncovered columns. DO NOT alter stars, primes or covers.
     *         Return to 3.
     *
     */
    // *******************************************************************************************
    void hungarian_method(const Eigen::MatrixXd &cost_matrix_arg, std::vector<Index> &optimal_assignments, const double _tol);

    //void reduce_cost(Eigen::MatrixXd &cost_matrix, double _infinity);

    //void find_zeros(const Eigen::MatrixXd &cost_matrix, Eigen::MatrixXi &zero_marks, double _tol);

    //bool check_assignment(const Eigen::MatrixXi &zero_marks, Eigen::VectorXi &col_covered);

    //int prime_zeros(const Eigen::MatrixXd &cost_matrix, Eigen::VectorXi &row_covered, Eigen::VectorXi &col_covered, Eigen::MatrixXi &zero_marks, double &min, Eigen::VectorXi &first_prime_zero);

    //int alternating_path(const Eigen::MatrixXd &cost_matrix, const Eigen::VectorXi &first_prime_zero, Eigen::MatrixXi &zero_marks, Eigen::VectorXi &row_covered, Eigen::VectorXi &col_covered);

    //int update_costs(const Eigen::VectorXi &row_covered, const Eigen::VectorXi &col_covered, const double min, Eigen::MatrixXd &cost_matrix);
  }

  //*******************************************************************************************
  ///Take a vector of doubles, and multiply by some factor that turns it into a vector of integers (within a tolerance)
  template <typename Derived>
  Eigen::Matrix<int,
        Derived::RowsAtCompileTime,
        Derived::ColsAtCompileTime>
  scale_to_int(const Eigen::MatrixBase<Derived> &val, double _tol = TOL) {

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
    if(valid_index(min_i))
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
    for(Index factor = 1; factor < 1000 && !within_tol; factor++) {
      Index i, j;
      tdubs = double(factor) * dubs;
      for(i = 0; i < dubs.rows(); i++) {
        for(j = 0; j < dubs.cols(); j++) {
          if(!almost_zero(round(tdubs(i, j)) - tdubs(i, j), _tol))
            break;
        }
        if(j < dubs.cols())
          break;
      }
      if(dubs.rows() <= i)
        within_tol = true;
    }

    if(within_tol) {
      for(Index i = 0; i < dubs.rows(); i++) {
        for(Index j = 0; j < dubs.cols(); j++) {
          ints(i, j) = round(tdubs(i, j));
        }
      }
    }

    return ints;
  }


}

namespace Eigen {
  template <typename Derived1, typename Derived2>
  inline
  bool almost_equal(const Eigen::MatrixBase<Derived1> &val1, const Eigen::MatrixBase<Derived2> &val2, double tol = CASM::TOL) {
    return CASM::almost_zero(val1 - val2, tol);
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
    return CASM::almost_zero(test_mat - test_mat.transpose(), test_tol);
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
    return CASM::almost_zero(test_mat - rev_mat.transpose(), test_tol);
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
