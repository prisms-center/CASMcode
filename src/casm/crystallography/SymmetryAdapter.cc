#include "casm/crystallography/SymmetryAdapter.hh"

namespace {
  /* SymOpMatrixType symop_to_matrix(const CASM::SymOp &op) { */
  /*   return op.matrix(); */
  /* } */
} // namespace

namespace CASM {
  namespace Adapter {

    SymOpType to_symop_type(const SymOpType &op) {
      return op;
    }

    SymOpMatrixType &get_matrix(SymOpType &op) {
      return op.matrix;
    }
    const SymOpMatrixType &get_matrix(const SymOpType &op) {
      return op.matrix;
    }

    SymOpTranslationType &get_translation(SymOpType &op) {
      return op.translation;
    }
    const SymOpTranslationType &get_translation(const SymOpType &op) {
      return op.translation;
    }

    SymOpTimeReversalType get_time_reversal(const SymOpType &op) {
      return op.is_time_reversal_active;
    }

  } // namespace Adapter
} // namespace CASM
