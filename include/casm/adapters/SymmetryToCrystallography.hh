#ifndef SYMMETRYTOCRYSTALLOGRAPHY_HH
#define SYMMETRYTOCRYSTALLOGRAPHY_HH

#include "casm/crystallography/SymmetryAdapter.hh"

namespace CASM {
  class SymOp;
  class SymGroup;
  namespace Adapter {
    SymOpType to_symop_type(const SymOp &op);
    /* SymGroupType to_symgroup_type(const SymGroup& group); */
  }
} // namespace CASM

#endif
