#ifndef CASM_enum_SuperConfigEnumInterface
#define CASM_enum_SuperConfigEnumInterface

#include "casm/app/enum/EnumInterface.hh"

namespace CASM {

  /// Interface for SuperConfigEnum
  class SuperConfigEnumInterface : public EnumInterfaceBase {
    CLONEABLE(SuperConfigEnumInterface)
  public:

    std::string desc() const override;

    std::string name() const override;

    void run(PrimClex &primclex, jsonParser const &json_options, jsonParser const &cli_options_as_json) const override;

  };

}

#endif
