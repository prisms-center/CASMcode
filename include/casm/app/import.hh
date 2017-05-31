#ifndef CASM_import
#define CASM_import

#include "casm/completer/Handlers.hh"
#include "casm/app/APICommand.hh"

namespace CASM {
  namespace Completer {

    class ImportOption : public OptionHandlerBase {

    public:

      using OptionHandlerBase::configtype;
      using OptionHandlerBase::configtype_opts;
      using OptionHandlerBase::settings_path;
      using OptionHandlerBase::input_str;


      ImportOption();

      const std::vector<fs::path> &pos_vec() const;

      const fs::path &batch_path() const;


    private:

      void initialize() override;

      std::vector<fs::path> m_pos_vec;

      fs::path m_batch_path;

    };
  }

  class ImportCommandImplBase;
  template<typename T> class ImportCommandImpl;

  /// 'casm import' command
  ///
  class ImportCommand :
    public APICommand<Completer::ImportOption> {

  public:

    using ImplBase = ImportCommandImplBase;
    template<typename ConfigType> using Impl = ImportCommandImpl<ConfigType>;

    static const std::string name;

    ImportCommand(const CommandArgs &_args, Completer::ImportOption &_opt);

    ~ImportCommand();

    int vm_count_check() const override;

    int help() const override;

    int desc() const override;

    int run() const override;

    // -- custom --

    ImportCommandImplBase &impl() const;

    void print_names(std::ostream &sout) const;

    jsonParser input() const;

  private:

    mutable std::unique_ptr<ImportCommandImplBase> m_impl;
  };

}

#endif
