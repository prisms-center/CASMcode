#ifndef UNITCELLCOORDTRAITS_HH
#define UNITCELLCOORDTRAITS_HH

#include "casm/clusterography/ElementSymApply.hh"
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

    typedef sym::CopyApplyWithPrim copy_apply_f_type;

    /* template <typename Base> */
    /* class CopyApplyType */
    /* { */
    /*     public: */
    /*     typedef CopyApplyWithPrim<Base> type; */
    /* }; */
  };
} // namespace CASM

#endif
