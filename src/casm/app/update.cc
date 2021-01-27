#include "casm/app/update.hh"

#include "casm/app/DBInterface.hh"
#include "casm/app/casm_functions.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/database/DatabaseTypes_impl.hh"
#include "casm/database/Import.hh"
#include "casm/database/Selection.hh"

// need to add specializations here
#include "casm/database/ConfigImport.hh"

namespace CASM {
namespace Completer {
UpdateOption::UpdateOption() : OptionHandlerBase("update") {}

void UpdateOption::initialize() {
  add_help_suboption();

  fs::path _default = "ALL";
  add_configlist_suboption(_default);

  add_configtype_suboption(traits<Configuration>::short_name,
                           DB::config_types_short());

  bool required = false;
  add_settings_suboption(required);
  add_input_suboption(required);

  return;
}

}  // namespace Completer

// -- class UpdateCommandImplBase --------------------------------------------

/// Defaults used if DataObject type doesn't matter or not given
class UpdateCommandImplBase {
 public:
  UpdateCommandImplBase(const UpdateCommand &cmd);

  virtual ~UpdateCommandImplBase() {}

  virtual int help() const;

  virtual int desc() const;

  virtual int run() const;

 protected:
  const UpdateCommand &m_cmd;
};

UpdateCommandImplBase::UpdateCommandImplBase(const UpdateCommand &cmd)
    : m_cmd(cmd) {}

int UpdateCommandImplBase::help() const {
  log() << std::endl << m_cmd.opt().desc() << std::endl;
  m_cmd.print_names(log());
  return 0;
}

int UpdateCommandImplBase::desc() const {
  help();

  log() << "Update calculation data for configurations specified by "
           "--selection. \n\n"

           "Calculated structures are mapped to the closest matching "
           "configuration \n"
           "consistent with the primitive crystal structure. \n\n"

           "For complete update options description for a particular config "
           "type, use\n"
           "--desc along with '--type <typename>'.\n\n";

  return 0;
}

int UpdateCommandImplBase::run() const {
  err_log() << "ERROR: No --type\n";
  m_cmd.print_names(err_log());
  return ERR_INVALID_ARG;
}

// -- template<typename DataObject> class UpdateCommandImpl -----------------

/// 'casm query' implementation, templated by type
///
/// This:
/// - holds a DB::InterfaceData object which stores dictionaries and selections
/// - provides the implementation for 'help' (i.e. print allowed import options)
/// - provides the implementation for 'run' (i.e. perform query)
///
template <typename DataObject>
class UpdateCommandImpl : public UpdateCommandImplBase {
 public:
  UpdateCommandImpl(const UpdateCommand &cmd) : UpdateCommandImplBase(cmd) {}

  int help() const override;

  int desc() const override;

  int run() const override;
};

template <typename DataObject>
int UpdateCommandImpl<DataObject>::help() const {
  log()
      << std::endl
      << m_cmd.opt().desc() << std::endl
      << "For complete options description for a particular config type, use:\n"
      << "  casm " << UpdateCommand::name << " --desc --type "
      << traits<DataObject>::short_name << "\n\n";
  return 0;
}

template <typename DataObject>
int UpdateCommandImpl<DataObject>::desc() const {
  log() << DB::Update<DataObject>::desc << std::endl;
  return 0;
}

template <typename DataObject>
int UpdateCommandImpl<DataObject>::run() const {
  return DB::Update<DataObject>::run(m_cmd.primclex(), m_cmd.input(),
                                     m_cmd.opt());
}

// -- class UpdateCommand ----------------------------------------------------

const std::string UpdateCommand::name = "update";

UpdateCommand::UpdateCommand(const CommandArgs &_args,
                             Completer::UpdateOption &_opt)
    : APICommand<Completer::UpdateOption>(_args, _opt) {}

UpdateCommand::~UpdateCommand() {}

int UpdateCommand::vm_count_check() const { return 0; }

int UpdateCommand::help() const { return impl().help(); }

int UpdateCommand::desc() const { return impl().desc(); }

int UpdateCommand::run() const { return impl().run(); }

UpdateCommandImplBase &UpdateCommand::impl() const {
  if (!m_impl) {
    if (vm().count("type")) {
      if (!opt().configtype_opts().count(opt().configtype())) {
        std::stringstream msg;
        msg << "--type " << opt().configtype() << " is not allowed for 'casm "
            << name << "'.";
        print_names(err_log());
        throw CASM::runtime_error(msg.str(), ERR_INVALID_ARG);
      }

      DB::for_config_type_short(
          opt().configtype(), DB::ConstructImpl<UpdateCommand>(m_impl, *this));
    } else {
      m_impl = notstd::make_unique<UpdateCommandImplBase>(*this);
    }
  }
  return *m_impl;
}

void UpdateCommand::print_names(std::ostream &sout) const {
  sout << "The allowed types are:\n";

  for (const auto &configtype : opt().configtype_opts()) {
    sout << "  " << configtype << std::endl;
  }
}

jsonParser UpdateCommand::input() const {
  if (count("settings")) {
    return jsonParser{opt().settings_path()};
  } else if (count("input")) {
    return jsonParser::parse(opt().input_str());
  }
  // use defaults
  return jsonParser();
}

}  // namespace CASM
