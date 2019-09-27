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

    template<typename SymGroupType>
    std::vector<SymOpMatrixType> symop_to_matrix(const SymGroupType &group) {
      std::vector<SymOpMatrixType> casted_group;
      casted_group.reserve(group.size());
      for(const auto &op : group) {
        casted_group.push_back(symop_to_matrix(op));
      }

      return casted_group;
    }
  } // namespace Adapter

} // namespace CASM

#endif
