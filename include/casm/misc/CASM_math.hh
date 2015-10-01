#ifndef CASM_MATH_HH
#define CASM_MATH_HH

#include "casm/CASM_global_definitions.hh"
#include "casm/container/Array.hh"
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
  inline
  bool almost_zero(const Eigen::MatrixXd &val, double tol = TOL) {
    return almost_zero(double(val.norm()), tol);
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
  //Return sign of number

  template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
  }

  // *******************************************************************************************

  // Works for signed and unsigned types
  template <typename IntType>
  IntType nchoosek(IntType n, IntType k) {
    assert(k <= n && 0 < k);
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

  double ran0(int &idum);

  // *******************************************************************************************
  /// Find greatest common factor
  int gcf(int i1, int i2);

  // *******************************************************************************************
  /// Find least common multiple
  int lcm(int i1, int i2);

  // *******************************************************************************************

  /* This creates Gaussians using the formula:
   * f(x) = a*e^(-(x-b)^2/(c^2))
   *************************************************************/

  double gaussian(double a, double x, double b, double c);

  // *******************************************************************************************


  Eigen::VectorXd eigen_vector_from_string(const std::string &tstr, const int &size);

  // *******************************************************************************************
  // calculates gcf(abs(i1),abs(i2)) and finds bezout coefficients p1 and p2 such that
  //             p1*i1 + p2*i2 = gcf(abs(i1),abs(i2));
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

    void reduce_cost(Eigen::MatrixXd &cost_matrix);

    void find_zeros(const Eigen::MatrixXd &cost_matrix, Eigen::MatrixXi &zero_marks, double _tol);

    bool check_assignment(const Eigen::MatrixXi &zero_marks, Eigen::VectorXi &col_covered);

    int prime_zeros(const Eigen::MatrixXd &cost_matrix, Eigen::VectorXi &row_covered, Eigen::VectorXi &col_covered, Eigen::MatrixXi &zero_marks, double &min, Eigen::VectorXi &first_prime_zero);

    int alternating_path(const Eigen::MatrixXd &cost_matrix, const Eigen::VectorXi &first_prime_zero, Eigen::MatrixXi &zero_marks, Eigen::VectorXi &row_covered, Eigen::VectorXi &col_covered);

    int update_costs(const Eigen::VectorXi &row_covered, const Eigen::VectorXi &col_covered, const double min, Eigen::MatrixXd &cost_matrix);
  }
}

namespace Eigen {
  template <typename Derived>
  inline
  bool almost_equal(const Eigen::MatrixBase<Derived> &val1, const Eigen::MatrixBase<Derived> &val2, double tol = CASM::TOL) {
    return CASM::almost_zero(double((val1 - val2).norm()), tol);
  }
}

#endif
