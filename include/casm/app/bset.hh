#ifndef CASM_bset
#define CASM_bset

#include "casm/app/APICommand.hh"

namespace CASM {
namespace Completer {
class BsetOption;
}

/// 'casm enum' implementation
class BsetCommand : public APICommand<Completer::BsetOption> {
 public:
  static const std::string name;

  BsetCommand(const CommandArgs &_args, Completer::BsetOption &_opt);

  ~BsetCommand();

  int vm_count_check() const override;

  int help() const override;

  int desc() const override;

  int run() const override;
};
}  // namespace CASM

#endif
