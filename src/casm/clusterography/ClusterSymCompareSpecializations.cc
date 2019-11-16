#include "casm/clusterography/ClusterSymCompare.hh"

namespace CASM {
  namespace xtal {
    class UnitCellCoord;
  }

  UnitCellCoord PrimPeriodicSymCompare<UnitCellCoord>::copy_apply_impl(SymOp const &op, UnitCellCoord obj) {
    return sym::copy_apply(op, obj, *m_prim);
  }

} // namespace CASM
