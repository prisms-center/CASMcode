#ifndef CASM_info_NeighborListInfoInterface
#define CASM_info_NeighborListInfoInterface

#include "casm/app/info/InfoInterface.hh"

namespace CASM {

/// Interface for neighbor list info
class NeighborListInfoInterface : public InfoInterfaceBase {
  CLONEABLE(NeighborListInfoInterface)
 public:
  std::string desc() const override;

  std::string name() const override;

  void run(jsonParser const &json_options,
           PrimClex const *primclex = nullptr) const override;
};

}  // namespace CASM

#endif
