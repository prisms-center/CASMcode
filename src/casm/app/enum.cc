#include <cstring>

#include "casm/app/enum.hh"
#include "casm/app/casm_functions.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/EnumeratorHandler.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/clex/SuperConfigEnum.hh"
#include "casm/completer/Handlers.hh"
#include "casm/casm_io/stream_io/container.hh"

namespace CASM {
  namespace Completer {

    EnumOption::EnumOption(): OptionHandlerBase("enum") {}

    void EnumOption::initialize() {
      bool required = false;

      m_desc.add_options()
      ("help,h", "Print help message.")
      ("desc",
       po::value<std::vector<std::string> >(&m_desc_vec)->multitoken()->zero_tokens()->value_name(ArgHandler::enummethod()),
       "Print extended usage description. "
       "Use '--desc MethodName [MethodName2...]' for detailed option description. "
       "Partial matches of method names will be included.")
      ("method,m", po::value<std::string>(&m_method)->value_name(ArgHandler::enummethod()), "Method to use: Can use number shortcuts in this option.")
      ("min", po::value<int>(&m_min_volume)->default_value(1), "Min volume")
      ("max", po::value<int>(&m_max_volume), "Max volume")
      ("filter",
       po::value<std::vector<std::string> >(&m_filter_strs)->multitoken()->value_name(ArgHandler::query()),
       "Filter configuration enumeration so that only configurations matching a "
       "'casm query'-type expression are recorded")
      ("all,a",
       po::bool_switch(&m_all_existing)->default_value(false),
       "Enumerate configurations for all existing supercells");

      add_verbosity_suboption();
      add_settings_suboption(required);
      add_input_suboption(required);
      add_scelnames_suboption();
      add_confignames_suboption();

      return;
    }
  }


  const std::string EnumCommand::name = "enum";

  EnumCommand::EnumCommand(const CommandArgs &_args, Completer::EnumOption &_opt) :
    APICommand<Completer::EnumOption>(_args, _opt),
    m_enumerators(nullptr) {}

  int EnumCommand::vm_count_check() const {
    if(!in_project()) {
      err_log().error("No casm project found");
      err_log() << std::endl;
      return ERR_NO_PROJ;
    }

    if(vm().count("method") != 1) {
      err_log() << "Error in 'casm enum'. The --method option is required." << std::endl;
      return ERR_INVALID_ARG;
    }

    if(vm().count("settings") + vm().count("input") == 2) {
      err_log() << "Error in 'casm enum'. The options --settings or --input may not both be chosen." << std::endl;
      return ERR_INVALID_ARG;
    }

    return 0;
  }

  int EnumCommand::help() const {
    log() << opt().desc() << std::endl;
    print_names(log(), enumerators());

    log() << "\nFor complete options description, use 'casm enum --desc MethodName'.\n\n";

    return 0;
  }

  int EnumCommand::desc() const {
    if(opt().desc_vec().size()) {
      log() << "\n";

      bool match = false;
      for(const auto &in_name : opt().desc_vec()) {
        for(const auto &e : enumerators()) {
          if(e.name().substr(0, in_name.size()) == in_name) {
            log() << e.help() << std::endl;
            match = true;
            break;
          }
        }
      }

      if(!match) {
        log() << "No match found. ";
        print_names(log(), enumerators());
      }

      return 0;
    }
    else {
      log() << "\n";
      log() << opt().desc() << std::endl;

      log() << "DESCRIPTION\n" << std::endl;

      log() << "  casm enum --settings input.json                                      \n"
            "  casm enum --input '{...JSON...}'                                     \n"
            "  - Input settings in JSON format to run an enumeration. The expected  \n"
            "    format is:                                                         \n"
            "\n"
            "    {\n"
            "      \"MethodName\": {\n"
            "        \"option1\" : ...,\n"
            "        \"option2\" : ...,\n"
            "         ...\n"
            "      }\n"
            "    }\n"
            "\n";

      print_names(log(), enumerators());

      log() << "\nFor complete options help for a particular method, \n"
            "use 'casm enum --desc MethodName'.\n\n";

      log() << "Custom enumerator plugins can be added by placing source code \n"
            "in the CASM project directory: \n"
            "  " << primclex().dir().enumerator_plugins() << " \n\n";

      return 0;
    }
  }

  int EnumCommand::run() const {
    jsonParser input;
    if(vm().count("settings")) {
      input = jsonParser {opt().settings_path()};
    }
    else if(vm().count("input")) {
      input = jsonParser::parse(opt().input_str());
    }

    auto lambda = [&](const EnumInterfaceBase & e) {
      return e.name().substr(0, opt().method().size()) == opt().method();
    };
    int count = std::count_if(enumerators().begin(), enumerators().end(), lambda);

    if(count == 1) {
      auto it = std::find_if(enumerators().begin(), enumerators().end(), lambda);
      return it->run(primclex(), input, opt());
    }
    else if(count < 1) {
      // allows for number aliasing
      try {
        int m = stoi(opt().method());
        if(m < std::distance(enumerators().begin(), enumerators().end()) &&
           m >= 0) {
          auto it = enumerators().begin();
          for(int k = 0; k < m; ++k) {
            ++it;
          }
          return it->run(primclex(), input, opt());
        }
      }
      catch(...) {}
      err_log() << "No match found for --method " << opt().method() << std::endl;
      print_names(err_log(), enumerators());
      return ERR_INVALID_ARG;
    }
    else if(count > 1) {
      err_log() << "Multiple matches found for --method " << opt().method() << std::endl;
      print_names(err_log(), enumerators());
      return ERR_INVALID_ARG;
    }
    throw std::runtime_error("Unknown error in EnumCommand::run");
  }

  const EnumeratorMap &EnumCommand::enumerators() const {
    if(!m_enumerators) {
      if(in_project()) {
        m_enumerators = &primclex().settings().enumerator_handler().map();
      }
      else {
        m_standard_enumerators = make_standard_enumerator_map();
        m_enumerators = m_standard_enumerators.get();
      }
    }
    return *m_enumerators;
  }

  void EnumCommand::print_names(std::ostream &sout, const EnumeratorMap &enumerators) const {
    sout << "The enumeration methods are:\n";
    int counter = 0;
    for(const auto &e : enumerators) {
      sout << "  " << counter << ") " << e.name() << std::endl;
      ++counter;
    }
  }

}

