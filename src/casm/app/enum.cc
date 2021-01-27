#include "casm/app/enum.hh"

#include <cstring>

#include "casm/app/DirectoryStructure.hh"
#include "casm/app/EnumeratorHandler.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/casm_functions.hh"
#include "casm/app/enum/EnumInterface.hh"
#include "casm/app/enum/standard_enumerator_interfaces.hh"
#include "casm/app/io/json_io_impl.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {
namespace Completer {

EnumOption::EnumOption() : EnumOptionBase("enum") { initialize(); }

void EnumOption::initialize() {
  bool required = false;

  // Standard options: --help, --desc, --settings, --input
  m_desc.add_options()(
      "help,h",
      "Print help message including a list of all available methods.")(
      "desc",
      po::value<std::vector<std::string> >(&m_desc_vec)
          ->multitoken()
          ->zero_tokens()
          ->value_name(ArgHandler::enummethod()),
      "Print extended usage description. "
      "Use '--desc [MethodName [MethodName2...]]' for detailed method "
      "descriptions. "
      "Partial matches of method names are acceptable.");
  add_settings_suboption(required);
  add_input_suboption(required);

  // `enum` method and method specific options:
  //    --method, --min, --max, --all, --scelnames, --confignames
  //
  // It is up to individual methods how to use these input, but prefer that CLI
  // options take
  //    precendence in collisions with JSON input provided by --settings or
  //    --input
  m_desc.add_options()(
      "method,m",
      po::value<std::string>(&m_method)->value_name(ArgHandler::enummethod()),
      "Method to use: Can use method name (including partial matches) or "
      "index.")(
      "min", po::value<int>(&m_min_volume),
      "Minimum volume supercell (integer, multiple of the prim volume)")(
      "max", po::value<int>(&m_max_volume),
      "Maximum volume supercell (integer, multiple of the prim volume)")(
      "all,a", po::bool_switch(&m_all_existing)->default_value(false),
      "Enumerate configurations for all existing supercells");
  add_scelnames_suboption();
  add_confignames_suboption();

  // Options that control enumeration in a generic way: --filter, --verbosity,
  // --dry-run
  m_desc.add_options()(
      "filter",
      po::value<std::string>(&m_filter_str)->value_name(ArgHandler::query()),
      "Filter configuration enumeration so that only configurations matching a "
      "'casm query'-type expression are recorded");
  add_verbosity_suboption();
  add_dry_run_suboption();

  return;
}
}  // namespace Completer

const std::string EnumCommand::name = "enum";

EnumCommand::EnumCommand(const CommandArgs &_args, Completer::EnumOption &_opt)
    : APICommand<Completer::EnumOption>(_args, _opt),
      m_enumerator_vector(nullptr) {}

int EnumCommand::vm_count_check() const {
  if (!in_project()) {
    err_log().error("No casm project found");
    err_log() << std::endl;
    return ERR_NO_PROJ;
  }

  if (vm().count("method") != 1) {
    err_log() << "Error in 'casm enum'. The --method option is required."
              << std::endl;
    return ERR_INVALID_ARG;
  }

  if (vm().count("settings") + vm().count("input") == 2) {
    err_log() << "Error in 'casm enum'. The options --settings or --input may "
                 "not both be chosen."
              << std::endl;
    return ERR_INVALID_ARG;
  }

  return 0;
}

int EnumCommand::help() const {
  log() << opt().desc() << std::endl;
  print_names(log(), enumerators());

  log() << "\nFor complete options description, use 'casm enum --desc "
           "MethodName'.\n\n";

  return 0;
}

int EnumCommand::desc() const {
  if (opt().desc_vec().size()) {
    log() << "\n";

    bool match = false;
    for (const auto &in_name : opt().desc_vec()) {
      for (const auto &interface_ptr : enumerators()) {
        if (interface_ptr->name().substr(0, in_name.size()) == in_name) {
          log() << interface_ptr->desc() << std::endl;
          match = true;
          break;
        }
      }
    }

    if (!match) {
      log() << "No match found. ";
      print_names(log(), enumerators());
    }

    return 0;
  } else {
    log() << "\n";
    log() << opt().desc() << std::endl;

    log() << "DESCRIPTION\n" << std::endl;

    log() << "  casm enum --settings input.json                                "
             "      \n"
             "  casm enum --input '{...JSON...}'                               "
             "      \n"
             "  - Input settings in JSON format to run an enumeration.         "
             "      \n\n";

    print_names(log(), enumerators());

    log() << "\nFor complete options help for a particular method, \n"
             "use 'casm enum --desc MethodName'.\n\n";

    log() << "Custom enumerator plugins can be added by placing source code \n"
             "in the CASM project directory: \n"
             "  "
          << primclex().dir().enumerator_plugins() << " \n\n";

    return 0;
  }
}

int EnumCommand::run() const {
  // Select and execute an enumeration method from this->enumerators() based on
  // --method value
  // - accepts partial name matches if only one partial match
  // - if no name match, attempt to interpret --method as index into
  // this->enumerators()
  // - otherwise provide a useful error message

  jsonParser json_options =
      make_json_input(opt());  // JSON from --input string or --settings file
  jsonParser cli_options_as_json{opt()};  // All CLI options as JSON object

  // find how many method names match the --method input value
  auto enumeration_method_name_matches =
      [&](notstd::cloneable_ptr<EnumInterfaceBase> const &interface_ptr) {
        return interface_ptr->name().substr(0, opt().method().size()) ==
               opt().method();
      };
  int count = std::count_if(enumerators().begin(), enumerators().end(),
                            enumeration_method_name_matches);

  if (count == 1) {
    auto it = std::find_if(enumerators().begin(), enumerators().end(),
                           enumeration_method_name_matches);
    (*it)->run(primclex(), json_options, cli_options_as_json);
    return 0;
  } else if (count < 1) {
    // Attempt to understand --method input as index into method list
    int method_index = -1;
    try {
      method_index = stoi(opt().method());
    } catch (...) {
      err_log() << "No match found for --method " << opt().method()
                << std::endl;
      print_names(err_log(), enumerators());
      return ERR_INVALID_ARG;
    }

    if (method_index < 0 || method_index >= enumerators().size()) {
      err_log() << "No match found for --method " << opt().method()
                << std::endl;
      print_names(err_log(), enumerators());
      return ERR_INVALID_ARG;
    }

    auto it = enumerators().begin();
    std::advance(it, method_index);
    (*it)->run(primclex(), json_options, cli_options_as_json);
    return 0;
  } else if (count > 1) {
    err_log() << "Multiple matches found for --method " << opt().method()
              << std::endl;
    print_names(err_log(), enumerators());
    return ERR_INVALID_ARG;
  }
  throw std::runtime_error("Unknown error in EnumCommand::run");
}

EnumInterfaceVector const &EnumCommand::enumerators() const {
  if (m_enumerator_vector == nullptr) {
    if (in_project()) {
      // include plugins and standard enumerator methods
      m_enumerator_vector = &primclex().settings().enumerator_handler().get();
    } else {
      // only show standard enumerator methods (can't run if not in a project,
      // but can see help messages)
      m_standard_enumerators = make_standard_enumerator_interfaces();
      m_enumerator_vector = &m_standard_enumerators;
    }
  }
  return *m_enumerator_vector;
}

void EnumCommand::print_names(std::ostream &sout,
                              EnumInterfaceVector const &enumerators) const {
  sout << "The enumeration methods are:\n";
  int counter = 0;
  for (auto const &interface_ptr : enumerators) {
    sout << "  " << counter << ") " << interface_ptr->name() << std::endl;
    ++counter;
  }
}

}  // namespace CASM
