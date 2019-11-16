#ifndef UNITCELLCOORDTRAITS_HH
#define UNITCELLCOORDTRAITS_HH

#include <string>

namespace CASM {

  template <typename T>
  struct traits;

  namespace xtal {
    class UnitCellCoord;
  }

  template <>
  struct traits<CASM::xtal::UnitCellCoord> {
    static const std::string name;
  };
} // namespace CASM

#endif
