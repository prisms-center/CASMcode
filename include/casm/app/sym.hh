#ifndef CASM_sym
#define CASM_sym

#include "casm/app/APICommand.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {
  /// 'casm enum' implementation
  class SymCommand : public APICommand<Completer::SymOption> {

  public:

    static const std::string name;

    SymCommand(const CommandArgs &_args, Completer::SymOption &_opt);

    int vm_count_check() const override;

    int help() const override;

    int desc() const override;

    int run() const override;

  };

}

#endif
