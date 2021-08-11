#include "casm/app/rm.hh"

#include "casm/app/DBInterface.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/database/DatabaseTypes.hh"
#include "casm/database/Import.hh"

// need to add specializations here
#include "casm/clex/RemoveSupercell.hh"
#include "casm/database/ConfigImport.hh"

namespace CASM {

namespace Completer {
RmOption::RmOption() : OptionHandlerBase("rm") {}

bool RmOption::force() const { return vm().count("force"); }

bool RmOption::data() const { return vm().count("data"); }

std::vector<std::string> const &RmOption::structure_properties() const {
  return m_structure_properties;
}

void RmOption::initialize() {
  add_help_suboption();
  add_names_suboption();
  add_selection_no_default_suboption();
  add_db_type_suboption(traits<Configuration>::short_name, DB::types_short());

  m_desc.add_options()(

      "data,d", "Remove calculation data only.")(

      "force,f",
      "Force remove including data and dependent objects (for --type=scel).")(

      "structure-properties",
      po::value<std::vector<std::string> >(&m_structure_properties)
          ->multitoken(),
      "Remove properties imported from specified structure files.\n");

  add_dry_run_suboption();

  return;
}
}  // namespace Completer

/*
  namespace DB {

    template<>
    class MakeInterfaceData<RmCommand> {

      template<typename DataObject>
      static InterfaceData<DataObject> eval(
        const RmCommand& cmd,
        InterfaceData<DataObject>& data) {

        // set 'selected' column
        if(cmd.in_project() && cmd.vm.count("selection")) {
          data.m_sel.resize(1);
          data.m_sel[0].reset(
            new Selection<DataObject>(
              cmd.primclex().db<DataObject>(),
              cmd.opt().selection_path()));
          data.m_ss << cmd.opt().selection_path();
        }
      }

    };

  }
*/

// -- RmCommandImplBase --------------------------------------------

/// Defaults used if DataObject type doesn't matter or not given
class RmCommandImplBase {
 public:
  RmCommandImplBase(const RmCommand &cmd);

  virtual ~RmCommandImplBase() {}

  virtual int help() const;

  virtual int desc() const;

  virtual int run() const;

 protected:
  const RmCommand &m_cmd;
};

RmCommandImplBase::RmCommandImplBase(const RmCommand &cmd) : m_cmd(cmd) {}

int RmCommandImplBase::help() const {
  log() << std::endl << m_cmd.opt().desc() << std::endl;
  m_cmd.print_names(log());
  m_cmd.print_config_names(log());
  return 0;
}

int RmCommandImplBase::desc() const {  // -- custom --
  help();

  log() << "Erase objects specified by --names or --selection, and data if "
           "applicable. \n\n"

           "Calculated structures are mapped to the closest matching "
           "configuration \n"
           "consistent with the primitive crystal structure. \n\n"

           "For complete update options description for a particular config "
           "type, use\n"
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
/// - provides the implementation for 'help' (i.e. print allowed import options)
/// - provides the implementation for 'run' (i.e. perform query)
///
template <typename DataObject>
class RmCommandImpl : public RmCommandImplBase {
 public:
  RmCommandImpl(const RmCommand &cmd) : RmCommandImplBase(cmd) {}

  int help() const override;

  int desc() const override;

  int run() const override;
};

template <typename DataObject>
int RmCommandImpl<DataObject>::help() const {
  log()
      << std::endl
      << m_cmd.opt().desc() << std::endl
      << "For complete options description for a particular config type, use:\n"
      << "  casm " << RmCommand::name << " --desc --type "
      << traits<DataObject>::short_name << "\n\n";
  return 0;
}

template <typename DataObject>
int RmCommandImpl<DataObject>::desc() const {
  log() << DB::Remove<DataObject>::desc() << std::endl;
  return 0;
}

template <typename DataObject>
int RmCommandImpl<DataObject>::run() const {
  return DB::Remove<DataObject>::run(m_cmd.primclex(), m_cmd.opt());
}

// -- class RmCommand ----------------------------------------------------

const std::string RmCommand::name = "rm";

RmCommand::RmCommand(const CommandArgs &_args, Completer::RmOption &_opt)
    : APICommand<Completer::RmOption>(_args, _opt) {}

RmCommand::~RmCommand() {}

int RmCommand::vm_count_check() const {
  if (!in_project()) {
    err_log().error("No casm project found");
    err_log() << std::endl;
    return ERR_NO_PROJ;
  }

  return 0;
}

int RmCommand::help() const { return impl().help(); }

int RmCommand::desc() const { return impl().desc(); }

int RmCommand::run() const { return impl().run(); }

RmCommandImplBase &RmCommand::impl() const {
  if (!m_impl) {
    if (in_project()) {
      if (DB::types_short().count(opt().db_type())) {
        DB::for_type_short(opt().db_type(),
                           DB::ConstructImpl<RmCommand>(m_impl, *this));
      } else {
        std::stringstream msg;
        msg << "--type " << opt().db_type() << " is not allowed for 'casm "
            << name << "'.";
        print_names(err_log());
        throw CASM::runtime_error(msg.str(), ERR_INVALID_ARG);
      }
    } else {
      m_impl = notstd::make_unique<RmCommandImplBase>(*this);
    }
  }
  return *m_impl;
}

void RmCommand::print_names(std::ostream &sout) const {
  sout << "The allowed types are:\n";

  for (const auto &db_type : opt().db_type_opts()) {
    sout << "  " << db_type << std::endl;
  }
}

void RmCommand::print_config_names(std::ostream &sout) const {
  sout << "The allowed types with --data option are:\n";

  for (const auto &db_type : opt().db_type_opts()) {
    if (DB::config_types_short().count(db_type)) {
      sout << "  " << db_type << std::endl;
    }
  }
}

}  // namespace CASM
