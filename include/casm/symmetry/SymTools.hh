#ifndef SYMTOOLS_HH
#define SYMTOOLS_HH

#include <vector>

namespace CASM {
  class SymGroup;
  class SymOp;
  namespace xtal {
    class Lattice;
    class UnitCellCoord;
    class Structure;
  }
  namespace sym {
    ///Returns the subgroup of the given group that keeps the lattice invariant
    SymGroup invariant_subgroup(const SymGroup &super_group, const xtal::Lattice &lat);

    template<typename OutputIt>
    OutputIt invariant_subgroup(const std::vector<SymOp> &super_group, const xtal::Lattice &lat, OutputIt result);

    /// \brief Apply SymOp to a UnitCellCoord
    xtal::UnitCellCoord &apply(const SymOp &op, xtal::UnitCellCoord &ucc, const xtal::Structure &prim);

    /// \brief Copy and apply SymOp to a UnitCellCoord
    xtal::UnitCellCoord copy_apply(const SymOp &op, const xtal::UnitCellCoord &ucc, const xtal::Structure &prim);
  }
}

#endif
