#include "casm/crystallography/SymmetryAdapter.hh"

namespace {
  /* SymOpMatrixType symop_to_matrix(const CASM::SymOp &op) { */
  /*   return op.matrix(); */
  /* } */
} // namespace

namespace CASM {
  namespace Adapter {
    SymOpMatrixType &get_matrix(SymOpType &op) {
      return std::get<0>(op);
    }
    const SymOpMatrixType &get_matrix(const SymOpType &op) {
      return std::get<0>(op);
    }

    SymOpTranslationType &get_translation(SymOpType &op) {
      return std::get<1>(op);
    }
    const SymOpTranslationType &get_translation(const SymOpType &op) {
      return std::get<1>(op);
    }

    SymOpTimeReversalType get_time_reversal(const SymOpType &op) {
      return std::get<2>(op);
    }

    SymOpType construct_sym_op(const SymOpMatrixType &sym_matrix,
                               const SymOpTranslationType &sym_translation,
                               SymOpTimeReversalType &sym_time_reversal) {
      return std::make_tuple(sym_matrix, sym_translation, sym_time_reversal);
    }

  } // namespace Adapter
} // namespace CASM
