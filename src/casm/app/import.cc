#include "casm/app/import.hh"
#include "casm/clex/Configuration.hh"
#include "casm/database/DatabaseTypes.hh"
#include "casm/database/Import.hh"
#include "casm/app/DBInterface.hh"
#include "casm/clex/PrimClex.hh"

// need to add specializations here
#include "casm/database/ConfigImport.hh"
#include "casm/database/DiffTransConfigImport.hh"


namespace CASM {

  namespace Completer {
    ImportOption::ImportOption(): OptionHandlerBase("import") {}

    const std::vector<fs::path> &ImportOption::pos_vec() const {
      return m_pos_vec;
    }

    const fs::path &ImportOption::batch_path() const {
      return m_batch_path;
    }

    void ImportOption::initialize() {
      add_help_suboption();

      m_desc.add_options()

      ("pos,p",
       po::value<std::vector<fs::path> >(&m_pos_vec)->multitoken()->value_name(ArgHandler::path()),
       "Path(s) to structure(s) being imported (multiple allowed, but no "
       "wild-card matching)")

      ("batch,b",
       po::value<fs::path>(&m_batch_path)->value_name(ArgHandler::path()),
       "Path to batch file, which should list one structure file path per line "
       "(can be used in combination with --pos)")

      ("data,d",
       "Attempt to extract calculation data from the enclosing "
       "directory of the structure files, if it is available")

      ("copy-additional-files",
       "Recursively copy other files from the same directory as the properties.calc.json file.");

      add_configtype_suboption(
        traits<Configuration>::short_name, DB::config_types_short());
      bool required = false;
      add_settings_suboption(required);
      add_input_suboption(required);


      return;
    }
  }


  // -- class ImportCommandImplBase --------------------------------------------

  /// Defaults used if DataObject type doesn't matter or not given
  class ImportCommandImplBase : public Logging {
  public:

    ImportCommandImplBase(const ImportCommand &cmd);

    virtual ~ImportCommandImplBase() {}

    virtual int help() const;

    virtual int desc() const;

    virtual int run() const;

  protected:

    const ImportCommand &m_cmd;
  };

  ImportCommandImplBase::ImportCommandImplBase(const ImportCommand &cmd) :
    Logging(cmd),
    m_cmd(cmd) {}

  int ImportCommandImplBase::help() const {
    log() << std::endl << m_cmd.opt().desc() << std::endl;
    m_cmd.print_names(log());
    return 0;
  }

  int ImportCommandImplBase::desc() const {  // -- custom --
    help();

    log() <<
          "Import structures or calculation data specified by --pos or --batch. \n\n"

          "Structures are mapped to the closest matching configuration consistent \n"
          "with the primitive crystal structure. If a JSON file is specified, it \n"
          "will be interpreted as a 'properties.calc.json' file.\n\n"

          "For complete import options description for a particular config type, use\n"
          "--desc along with '--type <typename>'.\n\n";

    return 0;
  }

  int ImportCommandImplBase::run() const {
    err_log() << "ERROR: No --type\n";
    m_cmd.print_names(err_log());
    return ERR_INVALID_ARG;
  }


  // -- template<typename DataObject> class ImportCommandImpl -----------------

  /// 'casm query' implementation, templated by type
  ///
  /// This:
  /// - provides the implementation for 'help' (i.e. print allowed import options)
  /// - provides the implementation for 'run' (i.e. perform query)
  ///
  template<typename DataObject>
  class ImportCommandImpl : public ImportCommandImplBase {
  public:
    ImportCommandImpl(const ImportCommand &cmd) :
      ImportCommandImplBase(cmd) {}

    int help() const override;

    int desc() const override;

    int run() const override;

  };

  template<typename DataObject>
  int ImportCommandImpl<DataObject>::help() const {
    log() << std::endl << m_cmd.opt().desc() << "\n\n"

          << "For complete options description for a particular config type, use:\n"
          "  casm " << ImportCommand::name << " --desc --type " << traits<DataObject>::short_name << "\n\n";
    return 0;
  }

  template<typename DataObject>
  int ImportCommandImpl<DataObject>::desc() const {
    log() << DB::Import<DataObject>::desc << std::endl;
    return 0;
  }

  template<typename DataObject>
  int ImportCommandImpl<DataObject>::run() const {
    return DB::Import<DataObject>::run(m_cmd.primclex(), m_cmd.input(), m_cmd.opt());
  }


  // -- class ImportCommand ----------------------------------------------------

  const std::string ImportCommand::name = "import";

  ImportCommand::ImportCommand(const CommandArgs &_args, Completer::ImportOption &_opt) :
    APICommand<Completer::ImportOption>(_args, _opt) {}

  ImportCommand::~ImportCommand() {}

  int ImportCommand::vm_count_check() const {
    if(!count("pos") && !count("batch")) {
      err_log() << "Error in 'casm import'. "
                "Use --pos or --batch to specify structures to be imported" << std::endl;
      return ERR_INVALID_ARG;
    }

    return 0;
  }

  int ImportCommand::help() const {
    return impl().help();
  }

  int ImportCommand::desc() const {
    return impl().desc();
  }

  int ImportCommand::run() const {
    return impl().run();
  }

  ImportCommandImplBase &ImportCommand::impl() const {
    if(!m_impl) {
      if(vm().count("type")) {
        if(!opt().configtype_opts().count(opt().configtype())) {
          std::stringstream msg;
          msg << "--type " << opt().configtype() << " is not allowed for 'casm " << name << "'.";
          print_names(err_log());
          throw CASM::runtime_error(msg.str(), ERR_INVALID_ARG);
        }

        DB::for_config_type_short(opt().configtype(), DB::ConstructImpl<ImportCommand>(m_impl, *this));
      }
      else {
        m_impl = notstd::make_unique<ImportCommandImplBase>(*this);
      }
    }
    return *m_impl;
  }

  void ImportCommand::print_names(std::ostream &sout) const {
    sout << "The allowed types are:\n";

    for(const auto &configtype : opt().configtype_opts()) {
      sout << "  " << configtype << std::endl;
    }
  }

  jsonParser ImportCommand::input() const {
    if(count("settings")) {
      return jsonParser {opt().settings_path()};
    }
    else if(count("input")) {
      return jsonParser::parse(opt().input_str());
    }
    // use defaults
    return jsonParser();
  }

}




