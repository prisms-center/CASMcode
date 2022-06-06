#ifndef CASM_enum_ConfigEnumInterfaceTemplate
#define CASM_enum_ConfigEnumInterfaceTemplate

#include "casm/app/enum/EnumInterface.hh"

namespace CASM {

/// Template for creating a ConfigEnum interface
class ConfigEnumInterfaceTemplate : public EnumInterfaceBase {
  CLONEABLE(ConfigEnumInterfaceTemplate)
 public:
  std::string desc() const override;

  std::string name() const override;

  void run(PrimClex &primclex, jsonParser const &json_options,
           jsonParser const &cli_options_as_json) const override;
};

}  // namespace CASM

#endif
