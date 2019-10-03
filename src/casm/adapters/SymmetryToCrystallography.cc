#include "casm/adapters/SymmetryToCrystallography.hh"
#include "casm/symmetry/SymOp.hh""

namespace CASM {
  namespace Adapter {
    SymOpType to_symop_type(const SymOp &op) {
      return SymOpType(op.matrix(), op.tau(), op.time_reversal());
    }
  }
}
