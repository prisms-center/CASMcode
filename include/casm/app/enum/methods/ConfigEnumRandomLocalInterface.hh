#ifndef CASM_enum_ConfigEnumRandomLocalInterface
#define CASM_enum_ConfigEnumRandomLocalInterface

#include "casm/app/enum/EnumInterface.hh"

namespace CASM {

  /// Interface for ConfigEnumRandomLocal
  class ConfigEnumRandomLocalInterface : public EnumInterfaceBase {
    CLONEABLE(ConfigEnumRandomLocalInterface)
  public:

    std::string desc() const override;

    std::string name() const override;

    void run(PrimClex &primclex, jsonParser const &json_options, jsonParser const &cli_options_as_json) const override;

  };

}

#endif
