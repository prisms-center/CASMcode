#include "casm/crystallography/SymType.hh"

namespace CASM {
  namespace xtal {
    const SymOpMatrixType &get_matrix(const SymOpType &op) {
      return op.matrix;
    }

    const SymOpTranslationType &get_translation(const SymOpType &op) {
      return op.translation;
    }

    SymOpTimeReversalType get_time_reversal(const SymOpType &op) {
      return op.is_time_reversal_active;
    }
  } // namespace xtal
} // namespace CASM
