#ifndef CASM_enum_ConfigEnumRandomOccupationsInterface
#define CASM_enum_ConfigEnumRandomOccupationsInterface

#include "casm/app/enum/EnumInterface.hh"

namespace CASM {

/// Interface for ConfigEnumRandomOccupations
class ConfigEnumRandomOccupationsInterface : public EnumInterfaceBase {
  CLONEABLE(ConfigEnumRandomOccupationsInterface)
 public:
  std::string desc() const override;

  std::string name() const override;

  void run(PrimClex &primclex, jsonParser const &json_options,
           jsonParser const &cli_options_as_json) const override;
};

}  // namespace CASM

#endif
