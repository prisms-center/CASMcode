#ifndef CASM_info_SupercellInfoInterface
#define CASM_info_SupercellInfoInterface

#include "casm/app/info/InfoInterface.hh"

namespace CASM {

/// Interface for supercell info
class SupercellInfoInterface : public InfoInterfaceBase {
  CLONEABLE(SupercellInfoInterface)
 public:
  std::string desc() const override;

  std::string name() const override;

  void run(jsonParser const &json_options, PrimClex const *primclex,
           fs::path root) const override;
};

}  // namespace CASM

#endif
