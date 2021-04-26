#include "casm/app/info.hh"

#include <cstring>

#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/casm_functions.hh"
#include "casm/app/info/InfoInterface.hh"
#include "casm/app/info/standard_info_method_interfaces.hh"
#include "casm/app/io/json_io_impl.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

namespace Completer {
InfoOption::InfoOption() : OptionHandlerBase("info") {}

void InfoOption::initialize() {
  bool required = false;
  add_settings_suboption(required);
  add_input_suboption(required);

  // Standard options: --help, --desc, --settings, --input
  m_desc.add_options()(

      "help,h",
      "Print help message including a list of all available methods.")(

      "desc",
      po::value<std::vector<std::string> >(&m_desc_vec)
          ->multitoken()
          ->zero_tokens()
          ->value_name(ArgHandler::infomethod()),
      "Print extended usage description. "
      "Use '--desc [MethodName [MethodName2...]]' for detailed method "
      "descriptions. Partial matches of method names are acceptable.")(

      "method,m",
      po::value<std::string>(&m_method)->value_name(ArgHandler::infomethod()),
      "Method to use: Can use method name (including partial matches) or "
      "index.");
}
}  // namespace Completer
}  // namespace CASM

namespace CASM {
// ///////////////////////////////////////
// 'sym' function for casm
//    (add an 'if-else' statement in casm.cpp to call this)

const std::string InfoCommand::name = "info";

InfoCommand::InfoCommand(const CommandArgs &_args, Completer::InfoOption &_opt)
    : APICommand<Completer::InfoOption>(_args, _opt),
      m_info_method_vector(nullptr) {}

int InfoCommand::vm_count_check() const {
  if (vm().count("method") != 1) {
    err_log() << "Error in 'casm info'. The --method option is required."
              << std::endl;
    return ERR_INVALID_ARG;
  }

  if (vm().count("settings") + vm().count("input") == 2) {
    err_log() << "Error in 'casm info'. The options --settings or --input may "
                 "not both be chosen."
              << std::endl;
    return ERR_INVALID_ARG;
  }

  return 0;
}

int InfoCommand::help() const {
  log() << opt().desc() << std::endl;
  print_names(log(), info_methods());

  log() << "\nFor complete options description, use 'casm info --desc "
           "MethodName'.\n\n";

  return 0;
}

int InfoCommand::desc() const {
  if (opt().desc_vec().size()) {
    log() << "\n";

    bool match = false;
    for (const auto &in_name : opt().desc_vec()) {
      for (const auto &interface_ptr : info_methods()) {
        if (interface_ptr->name().substr(0, in_name.size()) == in_name) {
          log() << interface_ptr->desc() << std::endl;
          match = true;
          break;
        }
      }
    }

    if (!match) {
      log() << "No match found. ";
      print_names(log(), info_methods());
    }

    return 0;
  } else {
    log() << "\n";
    log() << opt().desc() << std::endl;

    log() << "DESCRIPTION\n" << std::endl;

    log() << "  casm info --settings input.json                           \n"
             "  casm info --input '{...JSON...}'                          \n"
             "  - Input settings in JSON format.                          \n\n";

    print_names(log(), info_methods());

    log() << "\nFor complete options help for a particular method, \n"
             "use 'casm info --desc MethodName'.\n\n";

    // could provide help about plugins here

    return 0;
  }
}

int InfoCommand::run() const {
  // Select and execute an info method from this->info_methods() based on
  // --method value
  // - accepts partial name matches if only one partial match
  // - if no name match, attempt to interpret --method as index into
  // this->info_methods()
  // - otherwise provide a useful error message

  // JSON from --input string or --settings file
  jsonParser json_options = make_json_input(opt());

  // find how many method names match the --method input value
  auto method_name_matches =
      [&](notstd::cloneable_ptr<InfoInterfaceBase> const &interface_ptr) {
        return interface_ptr->name().substr(0, opt().method().size()) ==
               opt().method();
      };
  int count = std::count_if(info_methods().begin(), info_methods().end(),
                            method_name_matches);

  if (count == 1) {
    auto it = std::find_if(info_methods().begin(), info_methods().end(),
                           method_name_matches);
    (*it)->run(json_options, args().primclex, args().root);
    return 0;

  } else if (count < 1) {
    // Attempt to understand --method input as index into method list
    int method_index = -1;
    try {
      method_index = stoi(opt().method());
    } catch (...) {
      err_log() << "No match found for --method " << opt().method()
                << std::endl;
      print_names(err_log(), info_methods());
      return ERR_INVALID_ARG;
    }

    if (method_index < 0 || method_index >= info_methods().size()) {
      err_log() << "No match found for --method " << opt().method()
                << std::endl;
      print_names(err_log(), info_methods());
      return ERR_INVALID_ARG;
    }

    auto it = info_methods().begin();
    std::advance(it, method_index);

    (*it)->run(json_options, args().primclex, args().root);
    return 0;
  } else if (count > 1) {
    err_log() << "Multiple matches found for --method " << opt().method()
              << std::endl;
    print_names(err_log(), info_methods());
    return ERR_INVALID_ARG;
  }
  throw std::runtime_error("Unknown error in InfoCommand::run");
}

InfoInterfaceVector const &InfoCommand::info_methods() const {
  if (m_info_method_vector == nullptr) {
    // could add plugins here...

    m_standard_info_methods = make_standard_info_method_interfaces();
    m_info_method_vector = &m_standard_info_methods;
  }
  return *m_info_method_vector;
}

void InfoCommand::print_names(std::ostream &sout,
                              InfoInterfaceVector const &info_methods) const {
  sout << "The `casm info` methods are:\n";
  int counter = 0;
  for (auto const &interface_ptr : info_methods) {
    sout << "  " << counter << ") " << interface_ptr->name() << std::endl;
    ++counter;
  }
}

}  // namespace CASM
