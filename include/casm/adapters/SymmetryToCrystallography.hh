#ifndef SYMMETRYTOCRYSTALLOGRAPHY_HH
#define SYMMETRYTOCRYSTALLOGRAPHY_HH

#include "casm/crystallography/SymmetryAdapter.hh"

namespace CASM {
  class SymOp;
  class SymGroup;
  namespace xtal {
    namespace Adapter {
      template <>
      SymOpType to_symop_type<SymOp>(const SymOp &op);

      /* template <> */
      /* SymGroupType to_symgroup_type(const SymGroup& group); */

    } // namespace Adapter
  } // namespace xtal
} // namespace CASM

#endif
