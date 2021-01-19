#ifndef CASM_select
#define CASM_select

#include "casm/app/APICommand.hh"
#include "casm/app/DBInterface.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {
/*namespace Completer {

  /// Options set for `casm select`.
  class SelectOption : public OptionHandlerBase {

  public:

    using OptionHandlerBase::help_opt_vec;
    using OptionHandlerBase::selection_paths;
    using OptionHandlerBase::output_path;
    using OptionHandlerBase::db_type;
    using OptionHandlerBase::db_type_opts;

    SelectOption();

    const std::vector<std::string> &criteria_vec() const;

  private:

    void initialize() override;

    // vector necessary to allow --set/--set-on/--set-off with or without an
argument std::vector<std::string> m_criteria_vec;

  };

}*/
}

namespace CASM {

class SelectCommandImplBase;
template <typename T>
class SelectCommandImpl;

/// 'casm select' command
///
/// This class wraps a pointer to an implementation class
/// SelectCommandImpl<DataObject> which does most of the work. The
/// implementation classes are defined in the source file.
///
class SelectCommand : public APICommand<Completer::SelectOption> {
 public:
  using ImplBase = SelectCommandImplBase;
  template <typename ConfigType>
  using Impl = SelectCommandImpl<ConfigType>;

  static const std::string name;

  SelectCommand(const CommandArgs &_args, Completer::SelectOption &_opt);

  ~SelectCommand();

  int vm_count_check() const override;

  int help() const override;

  int desc() const override;

  int run() const override;

  // -- custom --

  SelectCommandImplBase &impl() const;

  void print_names(std::ostream &sout) const;

 private:
  mutable std::unique_ptr<SelectCommandImplBase> m_impl;
};
}  // namespace CASM

#endif
