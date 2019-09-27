#include "casm/crystallography/SymmetryAdapter.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/SymGroup.hh"

namespace CASM {
  namespace Adapter {
    SymOpMatrixType symop_to_matrix(const CASM::SymOp &op) {
      return op.matrix();
    }
  } // namespace Adapter
} // namespace CASM
