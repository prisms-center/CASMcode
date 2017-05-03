#include "casm/app/rm.hh"
#include "casm/app/DBInterface.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Configuration.hh"
#include "casm/database/Import.hh"
#include "casm/database/DatabaseTypeTraits.hh"

// need to add specializations here
#include "casm/clex/RemoveSupercell.hh"
#include "casm/database/ConfigImport.hh"

namespace CASM {

  namespace Completer {
    RmOption::RmOption(): OptionHandlerBase("rm") {}

    bool RmOption::force() const {
      return vm().count("force");
    }

    bool RmOption::dry_run() const {
      return vm().count("dry-run");
    }

    bool RmOption::data() const {
      return vm().count("data");
    }

    void RmOption::initialize() {
      add_help_suboption();
      add_names_suboption();
      add_selection_suboption();
      add_db_type_suboption(traits<Configuration>::short_name, DB::types_short());

      m_desc.add_options()
      ("dry-run,n", "Dry run")
      ("data,d", "Remove calculation data only.")
      ("force,f", "Force remove including data and dependent objects (for --type=scel).");

      return;
    }
  }


  // -- RmCommandImplBase --------------------------------------------

  /// Defaults used if DataObject type doesn't matter or not given
  class RmCommandImplBase : public Logging {
  public:

    RmCommandImplBase(const RmCommand &cmd);

    virtual ~RmCommandImplBase() {}

    virtual int help() const;

    virtual int desc() const;

    virtual int run() const;

  protected:

    const RmCommand &m_cmd;
  };


  RmCommandImplBase::RmCommandImplBase(const RmCommand &cmd) :
    Logging(cmd),
    m_cmd(cmd) {}

  int RmCommandImplBase::help() const {
    log() << std::endl << m_cmd.opt().desc() << std::endl;
    m_cmd.print_names(log());
    m_cmd.print_config_names(log());
    return 0;
  }

  int RmCommandImplBase::desc() const {  // -- custom --
    help();

    log() <<
          "Erase objects specified by --names or --selection, and data if applicable. \n\n"

          "Calculated structures are mapped to the closest matching configuration \n"
          "consistent with the primitive crystal structure. \n\n"

          "For complete update options description for a particular config type, use\n"
          "--desc along with '--type <typename>'.\n\n";

    return 0;
  }

  int RmCommandImplBase::run() const {
    err_log() << "ERROR: No --type\n";
    m_cmd.print_names(err_log());
    return ERR_INVALID_ARG;
  }


  // -- RmCommandImpl -----------------

  /// 'casm query' implementation, templated by type
  ///
  /// This:
  /// - holds a DB::InterfaceData object which stores dictionaries and selections
  /// - provides the implementation for 'help' (i.e. print allowed import options)
  /// - provides the implementation for 'run' (i.e. perform query)
  ///
  template<typename DataObject>
  class RmCommandImpl : public RmCommandImplBase {
  public:
    RmCommandImpl(const RmCommand &cmd) :
      RmCommandImplBase(cmd) {}

    int help() const override;

    int desc() const override;

    int run() const override;

  };

  template<typename DataObject>
  int RmCommandImpl<DataObject>::help() const {
    log() << std::endl << m_cmd.opt().desc() << std::endl
          << "For complete options description for a particular config type, use:\n"
          << "  casm " << RmCommand::name << " --desc --type " << traits<DataObject>::short_name << "\n\n";
    return 0;
  }

  template<typename DataObject>
  int RmCommandImpl<DataObject>::desc() const {
    log() << DB::Remove<DataObject>::desc << std::endl;
    return 0;
  }

  template<typename DataObject>
  int RmCommandImpl<DataObject>::run() const {
    return DB::Remove<DataObject>::run(m_cmd.primclex(), m_cmd.opt());
  }


  // -- class RmCommand ----------------------------------------------------

  RmCommand::RmCommand(const CommandArgs &_args, Completer::RmOption &_opt) :
    APICommand<Completer::RmOption>(_args, _opt) {}

  int RmCommand::vm_count_check() const {
    return 0;
  }

  int RmCommand::help() const {
    return impl().help();
  }

  int RmCommand::desc() const {
    return impl().desc();
  }

  int RmCommand::run() const {
    return impl().run();
  }

  RmCommandImplBase &RmCommand::impl() const {
    if(!m_impl) {
      if(vm().count("type")) {
        if(DB::types_short().count(opt().db_type())) {
          DB::for_type(opt().db_type(), DB::ConstructImpl<RmCommand>(m_impl, *this));
        }
        else {
          std::stringstream msg;
          msg << "--type " << opt().db_type() << " is not allowed for 'casm " << name << "'.";
          print_names(err_log());
          throw CASM::runtime_error(msg.str(), ERR_INVALID_ARG);
        }
      }
      else {
        m_impl = notstd::make_unique<RmCommandImplBase>(*this);
      }
    }
    return *m_impl;
  }

  void RmCommand::print_names(std::ostream &sout) const {
    sout << "The allowed types are:\n\n";

    for(const auto &db_type : opt().db_type_opts()) {
      sout << "  " << db_type << std::endl;
    }
  }

  void RmCommand::print_config_names(std::ostream &sout) const {
    sout << "The allowed types with --data option are:\n\n";

    for(const auto &db_type : opt().db_type_opts()) {
      if(DB::config_types_short().count(db_type)) {
        sout << "  " << db_type << std::endl;
      }
    }
  }

}


