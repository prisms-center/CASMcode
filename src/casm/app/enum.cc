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
      add_desc_vec_suboption();
      add_verbosity_suboption();

      // must have one
      add_settings_suboption(required);
      add_input_suboption(required);

      return;
    }
  }


  // ///////////////////////////////////////
  // 'enum' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int enum_command(const CommandArgs &args) {

    //casm enum --settings input.json
    //- enumerate supercells, configs, hop local configurations, etc.

    // Add standard enumerators here
    EnumeratorMap enumerators = make_enumerator_map();
    enumerators.insert(
      EnumInterface<ScelEnum>(),
      EnumInterface<ConfigEnumAllOccupations>(),
      EnumInterface<SuperConfigEnum>()
    );


    po::variables_map vm;

    Completer::EnumOption enum_opt;

    const fs::path &root = args.root;

    try {
      po::store(po::parse_command_line(args.argc, args.argv, enum_opt.desc()), vm); // can throw

      if(!vm.count("help") && !vm.count("desc")) {

        if(root.empty()) {
          args.err_log.error("No casm project found");
          args.err_log << std::endl;
          return ERR_NO_PROJ;
        }

        if(vm.count("settings") + vm.count("input") != 1) {
          std::cerr << "Error in 'casm enum'. One and only one of --settings or --input must be chosen." << std::endl;
          return ERR_INVALID_ARG;
        }
      }

      /** --help option
       */
      if(vm.count("help")) {
        std::cout << "\n";
        std::cout << enum_opt.desc() << std::endl;

        std::cout << "The current enumeration methods are:\n\n";

        for(const auto &e : enumerators) {
          std::cout << "  " << e.name() << std::endl;
        }

        std::cout << "\nFor complete options description, use 'casm enum --desc MethodName'.\n\n";

        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems

      if(vm.count("desc") && enum_opt.desc_vec().size()) {
        std::cout << "\n";

        bool match = false;
        for(const auto &in_name : enum_opt.desc_vec()) {
          for(const auto &e : enumerators) {
            if(e.name().substr(0, in_name.size()) == in_name) {
              std::cout << e.help() << std::endl;
              match = true;
            }
          }
        }

        if(!match) {
          std::cout << "No match found. The current enumeration methods are:\n\n";
          for(const auto &e : enumerators) {
            std::cout << "  " << e.name() << std::endl;
          }
        }

        return 0;
      }

      if(vm.count("desc")) {
        std::cout << "\n";
        std::cout << enum_opt.desc() << std::endl;

        std::cout << "DESCRIPTION\n" << std::endl;

        std::cout << "  casm enum --settings input.sjon                                      \n"
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
                  "\n"
                  "The current enumeration methods are:\n\n";

        for(const auto &e : enumerators) {
          std::cout << "  " << e.name() << std::endl;
        }

        std::cout << "\nFor complete options help, use 'casm enum --desc MethodName'.\n\n";

        return 0;
      }
    }
    catch(po::error &e) {
      std::cerr << enum_opt.desc() << std::endl;
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      std::cerr << enum_opt.desc() << std::endl;
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_UNKNOWN;

    }

    // If 'args.primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);
    const DirectoryStructure &dir = primclex.dir();
    const ProjectSettings &set = primclex.settings();

    jsonParser input;
    if(vm.count("settings")) {
      input = jsonParser {enum_opt.settings_path()};
    }
    else if(vm.count("input")) {
      input = jsonParser::parse(enum_opt.input_str());
    }
    if(input.size() != 1) {
      args.err_log.error("Expected a single attribute in the input");
      args.err_log << std::endl;
    }
    auto it = input.begin();
    int res = enumerators[it.name()].run(primclex, *it);

    std::cout << std::endl;

    return res;
  };

}

