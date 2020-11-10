#ifndef CASM_enum_ConfigEnumSiteDoFsInterface
#define CASM_enum_ConfigEnumSiteDoFsInterface

#include "casm/app/enum/EnumInterface.hh"

namespace CASM {

  /// Interface for ConfigEnumSiteDoFs
  class ConfigEnumSiteDoFsInterface : public EnumInterfaceBase {
    CLONEABLE(ConfigEnumSiteDoFsInterface)
  public:

    std::string desc() const override;

    std::string name() const override;

    void run(PrimClex &primclex, jsonParser const &json_options, jsonParser const &cli_options_as_json) const override;

  };

}

#endif
