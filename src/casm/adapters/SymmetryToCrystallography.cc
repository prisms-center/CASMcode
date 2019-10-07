#include "casm/crystallography/SymType.hh"
#include "casm/crystallography/SymTypeInterface.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {
  namespace xtal {
    /* template <> */
    /* SymOpMatrixType &get_matrix<SymOp>(SymOp &op) { */
    /*   return op.matrix(); */
    /* } */
    template <>
    const SymOpMatrixType &get_matrix<SymOp>(const SymOp &op) {
      return op.matrix();
    }

    /* template <> */
    /* SymOpTranslationType &get_translation<SymOp>(SymOp &op) { */
    /*   return op.tau(); */
    /* } */
    template <>
    const SymOpTranslationType &get_translation<SymOp>(const SymOp &op) {
      return op.tau();
    }

    template <>
    SymOpTimeReversalType get_time_reversal<SymOp>(const SymOp &op) {
      return op.time_reversal();
    }

  } // namespace xtal
} // namespace CASM
