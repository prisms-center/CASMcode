#include "casm/crystallography/SymmetryAdapter.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/SymGroup.hh"

namespace CASM {
  namespace Adapter {
    SymOpMatrixType symop_to_matrix(const CASM::SymOp &op) {
      return op.matrix();
    }
    std::vector<SymOpMatrixType> symop_to_matrix(const CASM::SymGroup &group) {
      std::vector<SymOpMatrixType> casted_group;
      casted_group.reserve(group.size());
      for(const auto &op : group) {
        casted_group.push_back(cast(op));
      }

      return casted_group;
    }
  } // namespace Adapter
} // namespace CASM
