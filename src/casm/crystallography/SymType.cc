#include "casm/crystallography/SymTypeInterface.hh"

namespace CASM {
  namespace xtal {

    template <>
    SymOpMatrixType &get_matrix<SymOpType>(SymOpType &op) {
      return op.matrix;
    }
    template <>
    const SymOpMatrixType &get_matrix<SymOpType>(const SymOpType &op) {
      return op.matrix;
    }

    template <>
    SymOpTranslationType &get_translation<SymOpType>(SymOpType &op) {
      return op.translation;
    }
    template <>
    const SymOpTranslationType &get_translation<SymOpType>(const SymOpType &op) {
      return op.translation;
    }

    template <>
    SymOpTimeReversalType get_time_reversal<SymOpType>(const SymOpType &op) {
      return op.is_time_reversal_active;
    }

    /* SymOpMatrixType &get_matrix(SymOpType &op) { */
    /*   return op.matrix; */
    /* } */
    /* const SymOpMatrixType &get_matrix(const SymOpType &op) { */
    /*   return op.matrix; */
    /* } */

    /* SymOpTranslationType &get_translation(SymOpType &op) { */
    /*   return op.translation; */
    /* } */
    /* const SymOpTranslationType &get_translation(const SymOpType &op) { */
    /*   return op.translation; */
    /* } */

    /* SymOpTimeReversalType get_time_reversal(const SymOpType &op) { */
    /*   return op.is_time_reversal_active; */
    /* } */

  } // namespace xtal
} // namespace CASM
