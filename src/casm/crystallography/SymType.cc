#include "casm/crystallography/SymType.hh"

namespace CASM {
  namespace xtal {
    const SymOpMatrixType &get_matrix(const SymOp &op) {
      return op.matrix;
    }

    const SymOpTranslationType &get_translation(const SymOp &op) {
      return op.translation;
    }

    SymOpTimeReversalType get_time_reversal(const SymOp &op) {
      return op.is_time_reversal_active;
    }
  } // namespace xtal
} // namespace CASM
