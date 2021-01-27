#ifndef CASM_global_eigen
#define CASM_global_eigen

#include "casm/external/Eigen/Dense"

namespace CASM {
typedef Eigen::MatrixXd::Index EigenIndex;
}

namespace Eigen {

typedef Matrix<long int, 3, 3> Matrix3l;
typedef Matrix<long int, 3, 1> Vector3l;
typedef Matrix<long int, Dynamic, Dynamic> MatrixXl;
typedef Matrix<long int, Dynamic, 1> VectorXl;

template <typename Derived>
std::istream &operator>>(std::istream &s, MatrixBase<Derived> &m) {
  for (int i = 0; i < m.rows(); ++i)
    for (int j = 0; j < m.cols(); j++) s >> m(i, j);
  return s;
}
}  // namespace Eigen

#endif
