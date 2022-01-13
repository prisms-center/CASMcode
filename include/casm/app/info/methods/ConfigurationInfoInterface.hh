#ifndef CASM_info_ConfigurationInfoInterface
#define CASM_info_ConfigurationInfoInterface

#include "casm/app/info/InfoInterface.hh"

namespace CASM {

/// Interface for configuration info
class ConfigurationInfoInterface : public InfoInterfaceBase {
  CLONEABLE(ConfigurationInfoInterface)
 public:
  std::string desc() const override;

  std::string name() const override;

  void run(jsonParser const &json_options, PrimClex const *primclex,
           fs::path root) const override;
};

}  // namespace CASM

#endif
