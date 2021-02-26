#ifndef CASM_symmetry_to_real
#define CASM_symmetry_to_real

#include "casm/container/multivector.hh"
#include "casm/external/Eigen/Core"
#include "casm/global/definitions.hh"

namespace CASM {

namespace SymRepTools_v2 {

template <typename T>
struct _RealType;

template <typename T>
using _Real = typename _RealType<T>::Type;

template <typename T>
struct _RealType<std::vector<T>> {
  using Type = std::vector<_Real<T>>;
};

template <typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
struct _RealType<Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>> {
  using Type = Eigen::Matrix<double, RowsAtCompileTime, ColsAtCompileTime>;
};

template <typename Derived>
_Real<Derived> to_real(Eigen::MatrixBase<Derived> const &mat) {
  return mat.real().template cast<double>();
}

template <typename T>
_Real<std::vector<T>> to_real(std::vector<T> const &vec) {
  std::vector<_Real<T>> result;
  result.reserve(vec.size());
  for (T const &el : vec) {
    result.emplace_back(to_real(el));
  }

  return result;
}

}  // namespace SymRepTools_v2

}  // namespace CASM

#endif
