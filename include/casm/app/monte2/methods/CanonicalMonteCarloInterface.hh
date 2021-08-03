#ifndef CASM_monte2_CanonicalMonteCarloInterface
#define CASM_monte2_CanonicalMonteCarloInterface

#include "casm/app/monte2/Monte2Interface.hh"

namespace CASM {

/// Interface for CanonicalMonteCarlo
class CanonicalMonteCarloInterface : public Monte2InterfaceBase {
  CLONEABLE(CanonicalMonteCarloInterface)
 public:
  std::string desc() const override;

  std::string name() const override;

  void run(PrimClex &primclex, jsonParser const &json_options,
           jsonParser const &cli_options_as_json) const override;
};

}  // namespace CASM

#endif
