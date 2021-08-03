#include "casm/app/monte2.hh"

#include <cstring>

#include "casm/app/DirectoryStructure.hh"
#include "casm/app/EnumeratorHandler.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/casm_functions.hh"
#include "casm/app/io/json_io_impl.hh"
#include "casm/app/monte2/Monte2Interface.hh"
#include "casm/app/monte2/standard_monte2_method_interfaces.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {
namespace Completer {

Monte2Option::Monte2Option() : OptionHandlerBase("monte2") { initialize(); }

void Monte2Option::initialize() {
  bool required = false;

  // Standard options: --help, --desc, --settings, --input
  m_desc.add_options()(
      "help,h",
      "Print help message including a list of all available methods.")(
      "desc",
      po::value<std::vector<std::string> >(&m_desc_vec)
          ->multitoken()
          ->zero_tokens()
          ->value_name(ArgHandler::montemethod()),
      "Print extended usage description. "
      "Use '--desc [MethodName [MethodName2...]]' for detailed method "
      "descriptions. "
      "Partial matches of method names are acceptable.");
  add_settings_suboption(required);
  add_input_suboption(required);

  // `monte2` method and method specific options:
  //    --method, --min, --max, --all, --scelnames, --confignames
  //
  // It is up to individual methods how to use these input, but prefer that CLI
  // options take
  //    precendence in collisions with JSON input provided by --settings or
  //    --input
  m_desc.add_options()(
      "method,m",
      po::value<std::string>(&m_method)->value_name(ArgHandler::montemethod()),
      "Method to use: Can use method name (including partial matches) or "
      "index.");

  // Options that control Monte Carlo calculations in a generic way: --verbosity
  add_verbosity_suboption();

  return;
}
}  // namespace Completer

const std::string Monte2Command::name = "monte2";

Monte2Command::Monte2Command(const CommandArgs &_args,
                             Completer::Monte2Option &_opt)
    : APICommand<Completer::Monte2Option>(_args, _opt),
      m_method_vector(nullptr) {}

int Monte2Command::vm_count_check() const {
  if (!in_project()) {
    err_log().error("No casm project found");
    err_log() << std::endl;
    return ERR_NO_PROJ;
  }

  if (vm().count("method") != 1) {
    err_log() << "Error in 'casm monte2'. The --method option is required."
              << std::endl;
    return ERR_INVALID_ARG;
  }

  if (vm().count("settings") + vm().count("input") == 2) {
    err_log()
        << "Error in 'casm monte2'. The options --settings or --input may "
           "not both be chosen."
        << std::endl;
    return ERR_INVALID_ARG;
  }

  return 0;
}

int Monte2Command::help() const {
  log() << opt().desc() << std::endl;
  print_names(log(), methods());

  log() << "\nFor complete options description, use 'casm monte2 --desc "
           "MethodName'.\n\n";

  return 0;
}

int Monte2Command::desc() const {
  if (opt().desc_vec().size()) {
    log() << "\n";

    bool match = false;
    for (const auto &in_name : opt().desc_vec()) {
      for (const auto &interface_ptr : methods()) {
        if (interface_ptr->name().substr(0, in_name.size()) == in_name) {
          log() << interface_ptr->desc() << std::endl;
          match = true;
          break;
        }
      }
    }

    if (!match) {
      log() << "No match found. ";
      print_names(log(), methods());
    }

    return 0;
  } else {
    log() << "\n";
    log() << opt().desc() << std::endl;

    log() << "DESCRIPTION\n" << std::endl;

    log() << "  casm monte2 --settings input.json\n"
             "  casm monte2 --input '{...JSON...}'\n"
             "  - Input settings in JSON format to run a Monte Carlo "
             "calculation.\n\n";

    print_names(log(), methods());

    log() << "\nFor complete options help for a particular method, \n"
             "use 'casm monte2 --desc MethodName'.\n\n";

    // log() << "Custom Monte Carlo method plugins can be added by placing \n"
    //          " source code in the CASM project directory: \n"
    //          "  "
    //       << primclex().dir().monte_method_plugins() << " \n\n";

    return 0;
  }
}

int Monte2Command::run() const {
  // Select and execute an Monte Carlo method from this->methods() based on
  // --method value
  // - accepts partial name matches if only one partial match
  // - if no name match, attempt to interpret --method as index into
  // this->methods()
  // - otherwise provide a useful error message

  jsonParser json_options =
      make_json_input(opt());  // JSON from --input string or --settings file
  jsonParser cli_options_as_json{opt()};  // All CLI options as JSON object

  // find how many method names match the --method input value
  auto monte_method_name_matches =
      [&](notstd::cloneable_ptr<Monte2InterfaceBase> const &interface_ptr) {
        return interface_ptr->name().substr(0, opt().method().size()) ==
               opt().method();
      };
  int count = std::count_if(methods().begin(), methods().end(),
                            monte_method_name_matches);

  if (count == 1) {
    auto it = std::find_if(methods().begin(), methods().end(),
                           monte_method_name_matches);
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
      print_names(err_log(), methods());
      return ERR_INVALID_ARG;
    }

    if (method_index < 0 || method_index >= methods().size()) {
      err_log() << "No match found for --method " << opt().method()
                << std::endl;
      print_names(err_log(), methods());
      return ERR_INVALID_ARG;
    }

    auto it = methods().begin();
    std::advance(it, method_index);
    (*it)->run(primclex(), json_options, cli_options_as_json);
    return 0;
  } else if (count > 1) {
    err_log() << "Multiple matches found for --method " << opt().method()
              << std::endl;
    print_names(err_log(), methods());
    return ERR_INVALID_ARG;
  }
  throw std::runtime_error("Unknown error in Monte2Command::run");
}

Monte2InterfaceVector const &Monte2Command::methods() const {
  if (m_method_vector == nullptr) {
    if (false /*in_project()*/) {
      // include plugins and standard methods
      // m_method_vector = &primclex().settings().monte_method_handler().get();
    } else {
      // only show standard method methods (can't run if not in a project,
      // but can see help messages)
      m_standard_methods = make_standard_monte2_method_interfaces();
      m_method_vector = &m_standard_methods;
    }
  }
  return *m_method_vector;
}

void Monte2Command::print_names(std::ostream &sout,
                                Monte2InterfaceVector const &methods) const {
  sout << "The Monte Carlo methods are:\n";
  int counter = 0;
  for (auto const &interface_ptr : methods) {
    sout << "  " << counter << ") " << interface_ptr->name() << std::endl;
    ++counter;
  }
}

}  // namespace CASM
