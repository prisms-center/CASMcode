#ifndef CASM_enum_ScelEnumInterface
#define CASM_enum_ScelEnumInterface

#include "casm/app/enum/EnumInterface.hh"

namespace CASM {

  /// Interface for ScelEnum
  class ScelEnumInterface : public EnumInterfaceBase {
    CLONEABLE(ScelEnumInterface)
  public:

    std::string desc() const override;

    std::string name() const override;

    void run(PrimClex &primclex, jsonParser const &json_options, jsonParser const &cli_options_as_json) const override;

  };

}

#endif
