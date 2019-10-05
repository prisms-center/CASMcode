#include "casm/adapters/SymmetryToCrystallography.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {
  namespace xtal {
    namespace Adapter {
      template <>
      xtal::Adapter::SymOpType to_symop_type<SymOp>(const SymOp &op) {
        return SymOpType(op.matrix(), op.tau(), op.time_reversal());
      }
    } // namespace Adapter
  } // namespace xtal
} // namespace CASM
