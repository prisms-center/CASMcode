#ifndef UNITCELLCOORDTRAITS_HH
#define UNITCELLCOORDTRAITS_HH

#include <string>

namespace CASM {
  template <typename Base>
  class CopyApplyWithPrim;

  template <typename T>
  struct traits;

  namespace xtal {
    class UnitCellCoord;
  }

  template <>
  struct traits<CASM::xtal::UnitCellCoord> {
    static const std::string name;
    template <typename Base>
    using CopyApplyType = CopyApplyWithPrim<Base>;

    /* template <typename Base> */
    /* class CopyApplyType */
    /* { */
    /*     public: */
    /*     typedef CopyApplyWithPrim<Base> type; */
    /* }; */
  };
} // namespace CASM

#endif
