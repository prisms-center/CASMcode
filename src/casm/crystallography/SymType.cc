#include "casm/crystallography/SymType.hh"

#include "casm/crystallography/Coordinate.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
namespace xtal {
SymOp operator*(const SymOp &LHS, const SymOp &RHS) {
  return SymOp(LHS.matrix * RHS.matrix,
               LHS.translation + LHS.matrix * RHS.translation,
               LHS.is_time_reversal_active !=
                   RHS.is_time_reversal_active);  // This is an XOR operation
}

//**********************************************************************************//

const SymOpMatrixType &get_matrix(const SymOp &op) { return op.matrix; }

const SymOpTranslationType &get_translation(const SymOp &op) {
  return op.translation;
}

SymOpTimeReversalType get_time_reversal(const SymOp &op) {
  return op.is_time_reversal_active;
}

}  // namespace xtal
}  // namespace CASM
