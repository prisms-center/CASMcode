#ifndef CASM_enum_ConfigEnumStrainInterface
#define CASM_enum_ConfigEnumStrainInterface

#include "casm/app/enum/EnumInterface.hh"

namespace CASM {

  /// Interface for ConfigEnumStrain
  class ConfigEnumStrainInterface : public EnumInterfaceBase {
    CLONEABLE(ConfigEnumStrainInterface)
  public:

    std::string desc() const override;

    std::string name() const override;

    void run(PrimClex &primclex, jsonParser const &json_options, jsonParser const &cli_options_as_json) const override;

  };

}

#endif
