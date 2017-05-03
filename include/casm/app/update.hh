#ifndef CASM_update
#define CASM_update

#include "casm/completer/Handlers.hh"
#include "casm/app/APICommand.hh"

namespace CASM {
  namespace Completer {

    class UpdateOption : public OptionHandlerBase {

    public:

      using OptionHandlerBase::configtype;
      using OptionHandlerBase::configtype_opts;
      using OptionHandlerBase::settings_path;
      using OptionHandlerBase::input_str;
      using OptionHandlerBase::selection_path;

      UpdateOption();

    private:

      void initialize() override;

    };
  }

  class UpdateCommandImplBase;
  template<typename T> class UpdateCommandImpl;

  /// 'casm update' command
  ///
  class UpdateCommand : public APICommand<Completer::UpdateOption> {

  public:

    using ImplBase = UpdateCommandImplBase;
    template<typename ConfigType> using Impl = UpdateCommandImpl<ConfigType>;

    static const std::string name;

    UpdateCommand(const CommandArgs &_args, Completer::UpdateOption &_opt);

    int vm_count_check() const override;

    int help() const override;

    int desc() const override;

    int run() const override;

    // -- custom --

    UpdateCommandImplBase &impl() const;

    void print_names(std::ostream &sout) const;

    jsonParser input() const;

  private:
    mutable std::unique_ptr<UpdateCommandImplBase> m_impl;
  };

}

#endif
