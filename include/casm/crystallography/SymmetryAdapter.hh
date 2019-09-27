#ifndef SYMMETRYADAPTER_HH
#define SYMMETRYADAPTER_HH

#include "casm/external/Eigen/Dense"
#include <vector>

namespace CASM {
  class SymOp;
  class SymGroup;
  namespace Adapter {
    typedef Eigen::Matrix3d SymOpMatrixType;

    SymOpMatrixType symop_to_matrix(const CASM::SymOp &op);
    std::vector<SymOpMatrixType> symop_to_matrix(const CASM::SymGroup &group);
    std::vector<SymOpMatrixType> symop_to_matrix(const std::vector<CASM::SymOp> &group);


  } // namespace Adapter

} // namespace CASM

#endif
