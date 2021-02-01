#include "casm/app/select.hh"

#include "casm/app/DBInterface_impl.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/clex/Configuration_impl.hh"
#include "casm/database/Database.hh"
#include "casm/database/DatabaseTypes.hh"
#include "casm/database/Selection.hh"

namespace CASM {

namespace Completer {
SelectOption::SelectOption() : OptionHandlerBase("select") {}

const std::vector<std::string> &SelectOption::criteria_vec() const {
  return m_criteria_vec;
}

void SelectOption::initialize() {
  add_general_help_suboption();
  add_selections_suboption();
  add_db_type_suboption(traits<Configuration>::short_name, DB::types_short());
  add_output_suboption("MASTER");

  m_desc.add_options()("json",
                       "Write JSON output (otherwise CSV, unless output "
                       "extension is '.json' or '.JSON')")(
      "subset",
      "Only write selected configurations to output. Can be used by itself or "
      "in conjunction with other options")(
      "xor", "Performs logical XOR on two configuration selections")(
      "not", "Performs logical NOT on configuration selection")(
      "or",
      "Write configurations selected in any of the input lists. Equivalent to "
      "logical OR")("and",
                    "Write configurations selected in all of the input lists. "
                    "Equivalent to logical AND")(
      "set-on",
      po::value<std::vector<std::string> >(&m_criteria_vec)
          ->multitoken()
          ->zero_tokens(),
      "Add configurations to selection if they meet specified criteria.  Call "
      "using 'casm select --set-on [\"criteria\"]'")(
      "set-off",
      po::value<std::vector<std::string> >(&m_criteria_vec)
          ->multitoken()
          ->zero_tokens(),
      "Remove configurations from selection if they meet specified criteria.  "
      "Call using 'casm select --set-off [\"criteria\"]'")(
      "set",
      po::value<std::vector<std::string> >(&m_criteria_vec)->multitoken(),
      "Create a selection of Configurations that meet specified criteria.  "
      "Call using 'casm select --set [\"criteria\"]'")("force,f",
                                                       "Overwrite output file");

  return;
}

}  // namespace Completer

// -- SelectCommandImplBase --------------------------------------------

/// Defaults used if DataObject type doesn't matter or not given
class SelectCommandImplBase {
 public:
  SelectCommandImplBase(const SelectCommand &cmd);

  virtual ~SelectCommandImplBase() {}

  virtual int help() const;

  virtual int desc() const;

  virtual int run() const;

 protected:
  const SelectCommand &m_cmd;
};

SelectCommandImplBase::SelectCommandImplBase(const SelectCommand &cmd)
    : m_cmd(cmd) {}

int SelectCommandImplBase::help() const {
  log() << std::endl << m_cmd.opt().desc() << std::endl;

  log()
      << "Use query commands to specify objects that should be selected or "
         "unselected.\n\n"

         "By default, the input and output selection is the MASTER selection, "
         "but one \n"
         "or more input selections may be specified via --selections (-c), and "
         "the output \n"
         "selection may be specified via --output (-o).\n\n"

         "For complete query options description, use '--help operators' or \n"
         "'--help properties' along with '--type <typename>'.\n\n"

         "The type of objects acted on is specified via --type (-t).\n\n";

  m_cmd.print_names(log());
  return 0;
}

int SelectCommandImplBase::desc() const { return help(); }

int SelectCommandImplBase::run() const {
  throw CASM::runtime_error("Unknown error in 'casm select'.", ERR_UNKNOWN);
}

// -- SelectCommandImpl -----------------

/// 'casm select' implementation, templated by type
///
/// This:
/// - holds a DB::InterfaceData object which stores dictionaries and selections
/// - provides the implementation for 'help' (i.e. print allowed query commands)
/// - provides the implementation for 'run' (i.e. --set, --set-on, --set-off,
/// --and, etc.)
///
template <typename DataObject>
class SelectCommandImpl : public SelectCommandImplBase {
 public:
  SelectCommandImpl(const SelectCommand &cmd);

  int help() const override;

  int desc() const override;

  int run() const override;

 private:
  int _count(std::string s) const { return m_cmd.vm().count(s); }

  const std::vector<fs::path> &_selection_paths() const {
    return m_cmd.opt().selection_paths();
  }

  fs::path _selection_paths(Index i) const { return _selection_paths()[i]; }

  const std::vector<std::string> &_criteria_vec() const {
    return m_cmd.opt().criteria_vec();
  }

  std::string _criteria() const;

  fs::path _output_path() const { return m_cmd.opt().output_path(); }

  void _set() const;

  void _set_on() const;

  void _set_off() const;

  void _and() const;

  void _or() const;

  void _xor() const;

  void _not() const;

  void _subset() const;

  void _write_input_stats() const;

  void _write_selection_stats(Index Ntot, const DB::Selection<DataObject> &sel,
                              Log &log, bool only_selected) const;

  const DataFormatterDictionary<DataObject> &_dict() const {
    return m_data.dict();
  }

  std::string _sel_str() const { return m_data.sel_str(); }

  double _sel_size() const { return m_data.sel_size(); }

  DB::Selection<DataObject> &_sel(Index i = 0) const { return m_data.sel(i); }

 private:
  // access dictionary and selections
  mutable DB::InterfaceData<DataObject> m_data;
  Index m_Ntot;
};

template <typename DataObject>
SelectCommandImpl<DataObject>::SelectCommandImpl(const SelectCommand &cmd)
    : SelectCommandImplBase(cmd),
      m_data(cmd),
      m_Ntot(m_cmd.primclex().template generic_db<DataObject>().size()) {}

template <typename DataObject>
int SelectCommandImpl<DataObject>::help() const {
  if (!m_cmd.opt().help_opt_vec().size()) {
    return SelectCommandImplBase::help();
  }

  for (const std::string &str : m_cmd.opt().help_opt_vec()) {
    if (str.empty()) {
      continue;
    }

    if (str[0] == 'o') {
      log() << "Available operators for use within selection criteria:"
            << std::endl;
      _dict().print_help(log(), DatumFormatterClass::Operator);
    } else if (str[0] == 'p') {
      log() << "Available property tags are currently:" << std::endl;
      _dict().print_help(log(), DatumFormatterClass::Property);
    }
    log() << std::endl;
  }
  log() << std::endl;
  return 0;
}

template <typename DataObject>
int SelectCommandImpl<DataObject>::desc() const {
  return help();
}

template <typename DataObject>
int SelectCommandImpl<DataObject>::run() const {
  _write_input_stats();

  log().begin_lap();
  if (_count("set")) {
    _set();
  } else if (_count("set-on")) {
    _set_on();
  } else if (_count("set-off")) {
    _set_off();
  } else if (_count("and")) {
    _and();
  } else if (_count("or")) {
    _or();
  } else if (_count("xor")) {
    _xor();
  } else if (_count("not")) {
    _not();
  } else {
    err_log() << "ERROR: No valid command recognized." << std::endl;
    help();
    return ERR_INVALID_ARG;
  }

  bool force_write = _count("force") || (_output_path().string() == "MASTER");
  if (fs::exists(_output_path()) && !force_write) {
    log() << "File " << _output_path()
          << " already exists. Use --force to force overwrite." << std::endl;
    return ERR_EXISTING_FILE;
  }

  bool only_selected = false;
  if (_count("subset")) {
    _subset();
    only_selected = true;
  }

  log() << "selection time: " << log().lap_time() << " (s)\n" << std::endl;

  log().write("Selection");
  _sel(0).write(_dict(), _output_path(), _count("json"), only_selected);

  log() << "write: " << _output_path() << "\n" << std::endl;

  log().custom("Output " + traits<DataObject>::short_name + " list",
               _output_path().string());
  _write_selection_stats(m_Ntot, _sel(0), log(), only_selected);

  log() << std::endl;

  return 0;
}

template <typename DataObject>
std::string SelectCommandImpl<DataObject>::_criteria() const {
  if (_criteria_vec().size() == 0) {
    return "";
  } else if (_criteria_vec().size() == 1) {
    return _criteria_vec()[0];
  } else {
    err_log()
        << "ERROR: Selection criteria must be a single string.  You provided "
        << _criteria_vec().size() << " strings:\n";
    for (const std::string &str : _criteria_vec())
      err_log() << "     - " << str << "\n";
    throw runtime_error("Invalid selection criteria", ERR_INVALID_ARG);
  }
}

template <typename DataObject>
void SelectCommandImpl<DataObject>::_set() const {
  log().custom("set", _criteria());
  _sel(0).set(_dict(), _criteria());
}

template <typename DataObject>
void SelectCommandImpl<DataObject>::_set_on() const {
  log().custom("set-on", _criteria());
  _sel(0).set(_dict(), _criteria(), true);
}

template <typename DataObject>
void SelectCommandImpl<DataObject>::_set_off() const {
  log().custom("set-off", _criteria());
  _sel(0).set(_dict(), _criteria(), false);
}

template <typename DataObject>
void SelectCommandImpl<DataObject>::_and() const {
  log().custom(std::string("and(") + _sel_str() + ")");
  for (int i = 1; i < _sel_size(); i++) {
    for (const auto &val : _sel(i).data()) {
      auto find_it = _sel(0).data().find(val.first);
      if (find_it != _sel(0).data().end()) {
        find_it->second = (find_it->second && val.second);
      } else {
        _sel(0).data()[val.first] = false;
      }
    }
  }
}

template <typename DataObject>
void SelectCommandImpl<DataObject>::_or() const {
  log().custom(std::string("or(") + _sel_str() + ")");
  for (int i = 1; i < _sel_size(); i++) {
    for (const auto &val : _sel(i).data()) {
      if (val.second) {
        _sel(0).data()[val.first] = true;
      }
    }
  }
}

template <typename DataObject>
void SelectCommandImpl<DataObject>::_xor() const {
  log().custom(_selection_paths(0).string() + " xor " +
               _selection_paths(1).string());
  for (const auto &val : _sel(1).data()) {
    // if not selected in second, use 'sel(0)' 'is_selected' value
    if (!val.second) {
      continue;
    }
    // else, if selected in second:

    // if not in 'sel(0)' insert selected
    auto find_it = _sel(0).data().find(val.first);
    if (find_it == _sel(0).data().end()) {
      _sel(0).data().insert(val);
    }
    // else, use opposite of sel(0) 'is_selected' value
    else {
      find_it->second = !find_it->second;
    }
  }
}

template <typename DataObject>
void SelectCommandImpl<DataObject>::_not() const {
  log().custom(std::string("not ") + _selection_paths(0).string());
  for (auto &value : _sel(0).data()) {
    value.second = !value.second;
  }
}

template <typename DataObject>
void SelectCommandImpl<DataObject>::_subset() const {
  auto it = _sel(0).data().cbegin();
  auto end = _sel(0).data().cend();
  while (it != end) {
    if (!it->second) {
      it = _sel(0).data().erase(it);
    } else {
      ++it;
    }
  }
}

template <typename DataObject>
void SelectCommandImpl<DataObject>::_write_input_stats() const {
  // ---- write starting stats ----
  log().custom("Input " + traits<DataObject>::short_name + " list",
               _selection_paths(0).string());
  _write_selection_stats(m_Ntot, _sel(0), log(), false);
  log() << std::endl;

  for (int i = 1; i < _selection_paths().size(); ++i) {
    log().custom("Input " + traits<DataObject>::short_name + " list",
                 _selection_paths(i).string());
    _write_selection_stats(m_Ntot, _sel(i), log(), false);
    log() << std::endl;
  }
}

template <typename DataObject>
void SelectCommandImpl<DataObject>::_write_selection_stats(
    Index Ntot, const DB::Selection<DataObject> &sel, Log &log,
    bool only_selected) const {
  auto Nselected = sel.selected_size();
  auto Ninclude = only_selected ? Nselected : sel.size();

  log << "# " + traits<DataObject>::short_name + "s in this project: " << Ntot
      << "\n";
  log << "# " + traits<DataObject>::short_name + "s included in this list: "
      << Ninclude << "\n";
  log << "# " + traits<DataObject>::short_name + "s selected in this list: "
      << Nselected << "\n";
}

// -- class SelectCommand ----------------------------------------------------

const std::string SelectCommand::name = "select";

SelectCommand::SelectCommand(const CommandArgs &_args,
                             Completer::SelectOption &_opt)
    : APICommand<Completer::SelectOption>(_args, _opt) {}

SelectCommand::~SelectCommand() {}

int SelectCommand::vm_count_check() const {
  if (!in_project()) {
    err_log().error("No casm project found");
    err_log() << std::endl;
    return ERR_NO_PROJ;
  }

  std::string cmd;
  std::vector<std::string> allowed_cmd = {"and",    "or",      "xor", "not",
                                          "set-on", "set-off", "set"};

  Index num_cmd(0);
  for (const std::string &cmd_str : allowed_cmd) {
    if (vm().count(cmd_str)) {
      num_cmd++;
      cmd = cmd_str;
    }
  }

  if (num_cmd > 1) {
    err_log() << "Error in 'casm select'. No more than one of the following "
                 "may be used: "
              << allowed_cmd << std::endl;
    return ERR_INVALID_ARG;
  } else if (vm().count("subset") && vm().count("selections") &&
             opt().selection_paths().size() != 1) {
    err_log()
        << "ERROR: 'casm select --subset' expects zero or one list as argument."
        << std::endl;
    return ERR_INVALID_ARG;
  }

  if (!vm().count("output") &&
      (cmd == "or" || cmd == "and" || cmd == "xor" || cmd == "not")) {
    err_log() << "ERROR: 'casm select --" << cmd
              << "' expects an --output file." << std::endl;
    return ERR_INVALID_ARG;
  }

  if ((cmd == "set-on" || cmd == "set-off" || cmd == "set") &&
      vm().count("selections") && opt().selection_paths().size() != 1) {
    err_log() << "Error in 'casm select " << cmd << "'. "
              << opt().selection_paths().size()
              << " selections were specified, "
                 "but no more than one selection is allowed (MASTER list is "
                 "used if no "
                 "other is specified)."
              << std::endl;
    return ERR_INVALID_ARG;
  }

  if (cmd == "xor" && opt().selection_paths().size() != 2) {
    err_log()
        << "ERROR: Option --xor requires exactly 2 selections as argument\n";
    return ERR_INVALID_ARG;
  }

  return 0;
}

int SelectCommand::help() const { return impl().help(); }

int SelectCommand::desc() const { return impl().desc(); }

int SelectCommand::run() const { return impl().run(); }

SelectCommandImplBase &SelectCommand::impl() const {
  if (!m_impl) {
    if (in_project()) {
      if (!opt().db_type_opts().count(opt().db_type())) {
        std::stringstream msg;
        msg << "--type " << opt().db_type() << " is not allowed for 'casm "
            << name << "'.";
        print_names(err_log());
        throw CASM::runtime_error(msg.str(), ERR_INVALID_ARG);
      }

      DB::for_type_short(opt().db_type(),
                         DB::ConstructImpl<SelectCommand>(m_impl, *this));
    } else {
      m_impl = notstd::make_unique<SelectCommandImplBase>(*this);
    }
  }
  return *m_impl;
}

void SelectCommand::print_names(std::ostream &sout) const {
  sout << "The allowed types are:\n";

  for (const auto &db_type : opt().db_type_opts()) {
    sout << "  " << db_type << std::endl;
  }
}

}  // namespace CASM
