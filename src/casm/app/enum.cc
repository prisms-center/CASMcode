#include <cstring>

#include "casm/app/casm_functions.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/clex/SuperConfigEnum.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

  namespace Completer {
    EnumOption::EnumOption(): OptionHandlerBase("enum") {}

    void EnumOption::initialize() {
      bool required = false;

      m_desc.add_options()
      ("help,h", "Print help message.")
      ("desc",
       po::value<std::vector<std::string> >(&m_desc_vec)->multitoken()->zero_tokens(),
       "Print extended usage description. "
       "Use '--desc MethodName [MethodName2...]' for detailed option description. "
       "Partial matches of method names will be included.")
      ("method,m", po::value<std::string>(&m_method), "Method to use: Can use number shortcuts in this option.")
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


  // ///////////////////////////////////////
  // 'enum' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int enum_command(const CommandArgs &args) {

    //casm enum --settings input.json
    //- enumerate supercells, configs, hop local configurations, etc.

    EnumeratorMap *enumerators;
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex *primclex;
    Completer::EnumOption enum_opt;
    po::variables_map &vm = enum_opt.vm();
    const fs::path &root = args.root;

    try {
      po::store(po::parse_command_line(args.argc, args.argv, enum_opt.desc()), vm); // can throw

      if(!vm.count("help") && !vm.count("desc")) {

        if(root.empty()) {
          args.err_log.error("No casm project found");
          args.err_log << std::endl;
          return ERR_NO_PROJ;
        }

        if(vm.count("method") != 1) {
          args.err_log << "Error in 'casm enum'. The --method option is required." << std::endl;
          return ERR_INVALID_ARG;
        }

        if(vm.count("settings") + vm.count("input") == 2) {
          args.err_log << "Error in 'casm enum'. The options --settings or --input may not both be chosen." << std::endl;
          return ERR_INVALID_ARG;
        }
      }

      if(!root.empty()) {
        primclex = &make_primclex_if_not(args, uniq_primclex);
        enumerators = &primclex->settings().enumerator_handler().map();
      }

      /** --help option
       */
      if(vm.count("help")) {
        args.log << "\n";
        args.log << enum_opt.desc() << std::endl;

        if(!root.empty()) {
          args.log << "The enumeration methods are:\n\n";

          int counter = 0;
          for(const auto &e : *enumerators) {
            args.log << "  " << counter << ") " << e.name() << std::endl;
            ++counter;
          }
        }

        args.log << "\nFor complete options description, use 'casm enum --desc MethodName'.\n\n";

        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems

      if(vm.count("desc") && enum_opt.desc_vec().size()) {
        args.log << "\n";

        bool match = false;
        for(const auto &in_name : enum_opt.desc_vec()) {
          for(const auto &e : *enumerators) {
            if(e.name().substr(0, in_name.size()) == in_name) {
              args.log << e.help() << std::endl;
              match = true;
            }
          }
        }

        if(!match) {

          if(!root.empty()) {
            args.log << "No match found. The enumeration methods are:\n\n";

            int counter = 0;
            for(const auto &e : *enumerators) {
              args.log << "  " << counter << ") " << e.name() << std::endl;
              ++counter;
            }
          }
        }

        return 0;
      }

      if(vm.count("desc")) {
        args.log << "\n";
        args.log << enum_opt.desc() << std::endl;

        args.log << "DESCRIPTION\n" << std::endl;

        args.log << "  casm enum --settings input.json                                      \n"
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

        if(!root.empty()) {
          args.log << "The enumeration methods are:\n\n";

          int counter = 0;
          for(const auto &e : *enumerators) {
            args.log << "  " << counter << ") " << e.name() << std::endl;
            ++counter;
          }


          args.log << "\nFor complete options help for a particular method, \n"
                   "use 'casm enum --desc MethodName'.\n\n";

          args.log << "Custom enumerator plugins can be added by placing source code \n"
                   "in the CASM project directory: \n"
                   "  " << primclex->dir().enumerator_plugins() << " \n\n";
        }

        return 0;
      }
    }
    catch(po::error &e) {
      args.err_log << enum_opt.desc() << std::endl;
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      args.err_log << enum_opt.desc() << std::endl;
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_UNKNOWN;
    }

    jsonParser input;
    if(vm.count("settings")) {
      input = jsonParser {enum_opt.settings_path()};
    }
    else if(vm.count("input")) {
      input = jsonParser::parse(enum_opt.input_str());
    }
    //int res = enumerators->find(enum_opt.method())->run(*primclex, input, enum_opt);
    auto lambda = [&](const EnumInterfaceBase & e) {
      return e.name().substr(0, enum_opt.method().size()) == enum_opt.method();
    };
    int count = std::count_if(enumerators->begin(), enumerators->end(), lambda);

    if(count == 1) {
      auto it = std::find_if(enumerators->begin(), enumerators->end(), lambda);

      return it->run(*primclex, input, enum_opt);
    }
    else if(count < 1) {
      // allows for number aliasing
      try {
        int m = stoi(enum_opt.method());
        if(m < std::distance(enumerators->begin(), enumerators->end()) &&
           m >= 0) {
          auto it = enumerators->begin();
          for(int k = 0; k < m; ++k) {
            ++it;
          }
          return it->run(*primclex, input, enum_opt);
        }
      }
      catch(...) {}
      args.log << "No match found. The enumeration methods are:\n\n";

      int counter = 0;
      for(const auto &e : *enumerators) {
        args.log << "  " << counter << ") " << e.name() << std::endl;
        ++counter;
      }
      return ERR_INVALID_ARG;

    }
    else if(count > 1) {
      args.err_log << "Multiple matches found for --method " << std::endl;
      return ERR_INVALID_ARG;
    }
    throw std::runtime_error("Unknown error in EnumCommand::run");

  };

}

