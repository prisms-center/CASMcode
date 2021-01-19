#ifndef CASM_MATH_HH
#define CASM_MATH_HH
#include "casm/global/definitions.hh"

/*
//Maybe we should transition to boost math library?
//#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>
*/

#include <cassert>
#include <cmath>
#include <complex>
#include <vector>

namespace CASM {

namespace CASM_TMP {

// --------------------

// Definitions for IfIntegralTol
template <typename tol_type, bool IsIntegral>
struct IfIntegralTol;

template <typename tol_type>
struct IfIntegralTol<tol_type, true> {
  IfIntegralTol(){};
  IfIntegralTol(tol_type){};
  tol_type tol() const { return 0; }
};

template <typename tol_type>
struct IfIntegralTol<tol_type, false> {
  IfIntegralTol(tol_type _tol) : m_tol(_tol){};
  tol_type tol() { return m_tol; }

 private:
  tol_type m_tol;
};

template <typename T>
using TypedTol = IfIntegralTol<T, std::is_integral<T>::value>;
// End of IfIntegralTol

// Definitions for MuchLessThan
template <typename value_type>
struct IntegralLessThan {
  IntegralLessThan(){};
  IntegralLessThan(value_type){};
  bool operator()(const value_type &A, const value_type &B) const {
    return A < B;
  };
};

template <typename value_type>
struct FloatingPointLessThan {
  FloatingPointLessThan(value_type _tol = TOL) : m_tol(_tol){};
  bool operator()(const value_type &A, const value_type &B) const {
    return A + m_tol < B;
  };

 private:
  value_type m_tol;
};

template <typename T>
using MuchLessThan =
    typename std::conditional<std::is_integral<T>::value, IntegralLessThan<T>,
                              FloatingPointLessThan<T> >::type;
// End of MuchLessThan

}  // namespace CASM_TMP

// *******************************************************************************************

// round function -- rounds to nearest whole number; half-integers are rounded
// away from zero
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
template <typename T, typename std::enable_if<std::is_floating_point<T>::value,
                                              T>::type * = nullptr>
inline bool almost_zero(const T &val, double tol = TOL) {
  return std::abs(val) < tol;
}

/// \brief If std::complex<T>, use std::abs(val) < tol;
template <typename T, typename std::enable_if<std::is_floating_point<T>::value,
                                              T>::type * = nullptr>
inline bool almost_zero(const std::complex<T> &val, double tol = TOL) {
  return std::abs(val) < tol;
}

/// \brief If T is integral, val == 0;
template <typename T, typename std::enable_if<std::is_integral<T>::value,
                                              T>::type * = nullptr>
inline bool almost_zero(const T &val, double tol = TOL) {
  return val == 0;
}

// *******************************************************************************************

/// \brief If T is not integral, use almost_zero(val1 - val2, tol);
template <typename T, typename std::enable_if<!std::is_integral<T>::value,
                                              T>::type * = nullptr>
bool almost_equal(const T &val1, const T &val2, double tol = TOL) {
  return almost_zero(val1 - val2, tol);
}

/// \brief If T is integral type, use val1 == val2;
template <typename T, typename std::enable_if<std::is_integral<T>::value,
                                              T>::type * = nullptr>
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
template <typename T, typename std::enable_if<std::is_floating_point<T>::value,
                                              T>::type * = nullptr>
bool compare(const T &A, const T &B, double tol) {
  if (!almost_equal(A, B, tol)) {
    return A < B;
  }
  return false;
}

struct FloatCompare {
  double tol;

  FloatCompare(double _tol) : tol(_tol) {}

  template <typename T>
  bool operator()(const T &A, const T &B) const {
    return compare(A, B, tol);
  }
};

/// \brief Floating point lexicographical comparison with tol
template <class InputIt1, class InputIt2>
bool float_lexicographical_compare(InputIt1 first1, InputIt1 last1,
                                   InputIt2 first2, InputIt2 last2,
                                   double tol) {
  FloatCompare compare(tol);
  return std::lexicographical_compare(first1, last1, first2, last2, compare);
}

// *******************************************************************************************
// Return sign of integral number
template <typename T, typename std::enable_if<std::is_integral<T>::value,
                                              T>::type * = nullptr>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

// *******************************************************************************************
// Return sign of floating-point number
template <typename T, typename std::enable_if<std::is_floating_point<T>::value,
                                              T>::type * = nullptr>
int float_sgn(T val, double compare_tol = TOL) {
  T zeroval(0);
  if (compare(zeroval, val, compare_tol)) {
    return 1;
  }

  else if (compare(val, zeroval, compare_tol)) {
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
  if (n < k) return 0;

  assert(0 <= k);
  if (n < 2 * k) k = n - k;

  n -= k;

  IntType result(1);
  for (IntType i = 1; i < k + 1; i++) {
    result *= n + i;
    result /= i;
  }

  return result;
}

// *******************************************************************************************

// Works for signed and unsigned types
template <typename IntType>
std::vector<IntType> index_to_kcombination(IntType ind, IntType k) {
  std::vector<IntType> result;
  result.reserve(k);
  IntType n;
  IntType big, bigger;
  for (; k > 0; --k) {
    n = k;
    bigger = 1;
    while (bigger <= ind) {
      big = bigger;
      bigger = nchoosek(++n, k);
    }
    result.push_back(n - 1);
    ind -= big;
  }
  return result;
}

// *******************************************************************************************
/// \brief Computes the Damerescau-Levenshtein distance
///   -- the number of edits (deletions, insertions, transpositions) to go from
///   string 'a' to string 'b'
int dl_string_dist(const std::string &a, const std::string &b);

// *******************************************************************************************

double ran0(int &idum);

// *******************************************************************************************
/// Find greatest common factor
int gcf(int i1, int i2);

// *******************************************************************************************
/// Find least common multiple
int lcm(int i1, int i2);

// *******************************************************************************************
/// \brief Calculate greatest common factor of two integers, and bezout
/// coefficients
///
/// \returns greatest common factor
////
/// \param i1,i2 two integers for which to find greatest common factor
/// \param[out] p1, p2 bezout coefficients such that p1*i1 + p2*i2 =
/// gcf(abs(i1),abs(i2));
///
/// \see smith_normal_form
///
template <typename IntType>
IntType extended_gcf(IntType i1, IntType i2, IntType &p1, IntType &p2) {
  using std::swap;

  IntType s1 = sgn(i1);
  IntType s2 = sgn(i2);

  IntType a = i1, b = i2;
  p1 = 1;
  p2 = 0;
  IntType tp1(0), tp2(1), quotient;

  i1 = std::abs(i1);
  i2 = std::abs(i2);

  while (i2 != 0) {
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

// evaluates Gaussians using the formula:
// f(x) = a*e^(-(x-b)^2/(c^2))
double gaussian(double a, double x, double b, double c);

// calculates Gaussian moments given by the integral:
// m = \int_{\infty}^{\infty} dx
// x^pow*exp[-x^2/(2*sigma^2)]/(\sqrt(2*\pi)*sigma)
double gaussian_moment(int expon, double sigma);

// calculates Gaussian moments given by the integral:
// m = \int_{\infty}^{\infty} dx
// x^pow*exp[-(x-x0)^2/(2*sigma^2)]/(\sqrt(2*\pi)*sigma)
double gaussian_moment(int expon, double sigma, double x0);

// *******************************************************************************************
// finds rational number that approximates 'val' to within 'tol' --
//      --  almost_int(double(denominator)*val, tol) == true OR
//      almost_int(val/double(numerator), tol) == true OR
//      --  denominator is always positive
//      --  sign(numerator)=sign(val)
//      --  if(almost_zero(val,tol)) --> numerator = 0, denominator = 1
// *******************************************************************************************

void nearest_rational_number(double val, long &numerator, long &denominator,
                             double tol = TOL);

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
/*
 * Returns string representation of an integer that has the same number of
 * characters as its maximum allowed value. To fix the length of the resulting
 * string, the character 'prepend_char' is prepended. When prepend_char is '0',
 * the resulting strings have sequential numerical values when ordered
 * alphabetically.
 */
std::string to_sequential_string(Index i, Index max_i, char prepend_char = '0');

// *******************************************************************************************
// John G 010413
int mod(int a, int b);

// *******************************************************************************************

double cuberoot(double number);

}  // namespace CASM

#endif
