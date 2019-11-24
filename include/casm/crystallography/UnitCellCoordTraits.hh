#ifndef UNITCELLCOORDTRAITS_HH
#define UNITCELLCOORDTRAITS_HH

#include "casm/symmetry/ElementSymApply.hh"
#include <string>

namespace CASM {
  template <typename Base>
  class CopyApplyWithPrim_crtp;

  template <typename T>
  struct traits;

  namespace xtal {
    class UnitCellCoord;
  }

  template <>
  struct traits<CASM::xtal::UnitCellCoord> {
    static const std::string name;
    template <typename Base>
    using CopyApplyType = CopyApplyWithPrim_crtp<Base>;

    typedef sym::CopyApplyWithPrim_f copy_apply_f_type;
  };
} // namespace CASM

#endif
