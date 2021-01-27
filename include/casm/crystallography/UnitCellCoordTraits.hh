#ifndef UNITCELLCOORDTRAITS_HH
#define UNITCELLCOORDTRAITS_HH

#include <string>

namespace CASM {
namespace sym {
class CopyApplyWithPrim_f;
}

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
  using copy_apply_crtp_type = CopyApplyWithPrim_crtp<Base>;

  typedef sym::CopyApplyWithPrim_f copy_apply_f_type;
};
}  // namespace CASM

#endif
