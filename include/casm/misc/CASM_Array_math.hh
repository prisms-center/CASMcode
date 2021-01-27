#ifndef CASM_Array_math
#define CASM_Array_math

#include "casm/global/definitions.hh"
//#include "casm/global/eigen.hh"
#include "casm/container/Array.hh"
#include "casm/misc/CASM_Eigen_math.hh"

/*
//Maybe we should transition to boost math library?
//#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>
*/

namespace CASM {

// //////////////////////////////////////////
// //////////////////////////////////////////
// Array Function declarations:

// *******************************************************************************************

// returns 'i' if 'input' is equivalent to 'unique[i]', w.r.t. permutation of
// the equivalent elements of 'input'. equivalent elements are specified by
// 'ind_equiv' if 'input' specifies a new combination of integers, unique.size()
// is returned
Index which_unique_combination(const Array<Index> &input,
                               const Array<Index>::X2 &unique,
                               const Array<Index>::X2 &ind_equiv);

// *******************************************************************************************

Index which_unique_combination(const Array<Index> &input,
                               const Array<Index>::X2 &unique);

// *******************************************************************************************
/// Find least common multiple
int lcm(const Array<int> &series);

// *******************************************************************************************

ReturnArray<Array<int> > get_prime_factors(int target);

// *******************************************************************************************

template <typename IntType>
IntType multinomial_coeff(const Array<IntType> &exponents) {
  IntType tcoeff(1), tsum(0);
  for (Index i = 0; i < exponents.size(); i++) {
    tsum += exponents[i];
    tcoeff *= nchoosek(tsum, exponents[i]);
  }
  return tcoeff;
}

// *******************************************************************************************
// get multinomial coefficient for only a subset of the coeffs
template <typename IntType>
IntType multinomial_coeff(const Array<IntType> &exponents,
                          const Array<Index> &sublist) {
  IntType tcoeff(1), tsum(0);
  for (Index i = 0; i < sublist.size(); i++) {
    tsum += exponents[sublist[i]];
    tcoeff *= nchoosek(tsum, exponents[sublist[i]]);
  }
  return tcoeff;
}
// ************************************************************
template <typename T>
bool almost_equal(const Array<T> &A, const Array<T> &B, double tol = TOL) {
  if (A.size() != B.size()) return false;
  for (Index i = 0; i < A.size(); i++) {
    if (!almost_equal(A[i], B[i], tol)) return false;
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
  for (Index i = 0; i < arr.size(); i++) {
    result.push_back(result[i] + arr[i]);
  }
  return result;
}

// ************************************************************

template <typename Derived>
ReturnArray<Index> partition_distinct_values(
    const Eigen::MatrixBase<Derived> &value, double tol = TOL) {
  Array<Index> subspace_dims;
  Index last_i = 0;
  for (Index i = 1; i < value.size(); i++) {
    if (!almost_equal(value[last_i], value[i], tol)) {
      subspace_dims.push_back(i - last_i);
      last_i = i;
    }
  }
  subspace_dims.push_back(value.size() - last_i);
  return subspace_dims;
}

}  // namespace CASM

#endif
