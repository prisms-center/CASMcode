#ifndef CASM_info_PrimInfoInterface
#define CASM_info_PrimInfoInterface

#include "casm/app/info/InfoInterface.hh"

namespace CASM {

/// Interface for prim info
class PrimInfoInterface : public InfoInterfaceBase {
  CLONEABLE(PrimInfoInterface)
 public:
  std::string desc() const override;

  std::string name() const override;

  void run(jsonParser const &json_options, PrimClex const *primclex,
           fs::path root) const override;
};

}  // namespace CASM

#endif
