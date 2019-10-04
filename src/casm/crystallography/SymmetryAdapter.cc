#include "casm/crystallography/SymmetryAdapter.hh"

namespace {
  /* SymOpMatrixType symop_to_matrix(const CASM::SymOp &op) { */
  /*   return op.matrix(); */
  /* } */
} // namespace

namespace CASM {
  namespace Adapter {

    template<>
    SymOpType to_symop_type<SymOpType>(const SymOpType &op) {
      return op;
    }

    template<>
    SymOpMatrixType &get_matrix<SymOpType>(SymOpType &op) {
      return op.matrix;
    }
    template<>
    const SymOpMatrixType &get_matrix<SymOpType>(const SymOpType &op) {
      return op.matrix;
    }

    template<>
    SymOpTranslationType &get_translation<SymOpType>(SymOpType &op) {
      return op.translation;
    }
    template<>
    const SymOpTranslationType &get_translation<SymOpType>(const SymOpType &op) {
      return op.translation;
    }

    template<>
    SymOpTimeReversalType get_time_reversal<SymOpType>(const SymOpType &op) {
      return op.is_time_reversal_active;
    }

  } // namespace Adapter
} // namespace CASM
