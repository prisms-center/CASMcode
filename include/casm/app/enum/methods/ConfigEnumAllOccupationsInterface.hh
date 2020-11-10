#ifndef CASM_enum_ConfigEnumAllOccupationsInterface
#define CASM_enum_ConfigEnumAllOccupationsInterface

#include "casm/app/enum/EnumInterface.hh"

namespace CASM {

  /// Interface for ConfigEnumAllOccupations
  class ConfigEnumAllOccupationsInterface : public EnumInterfaceBase {
    CLONEABLE(ConfigEnumAllOccupationsInterface)
  public:

    std::string desc() const override;

    std::string name() const override;

    void run(PrimClex &primclex, jsonParser const &json_options, jsonParser const &cli_options_as_json) const override;

  };

}

#endif
