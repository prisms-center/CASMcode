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

    EnumeratorMap *enumerators;
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex *primclex;
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
          args.err_log << "Error in 'casm enum'. One and only one of --settings or --input must be chosen." << std::endl;
          return ERR_INVALID_ARG;
        }
      }

      if(!root.empty()) {
        primclex = &make_primclex_if_not(args, uniq_primclex);
        enumerators = &primclex->enumerator_handler().map();
      }

      /** --help option
       */
      if(vm.count("help")) {
        args.log << "\n";
        args.log << enum_opt.desc() << std::endl;

        if(!root.empty()) {
          args.log << "The enumeration methods are:\n\n";
        
          for(const auto &e : *enumerators) {
            args.log << "  " << e.name() << std::endl;
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
          
            for(const auto &e : *enumerators) {
              args.log << "  " << e.name() << std::endl;
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

          for(const auto &e : *enumerators) {
            args.log << "  " << e.name() << std::endl;
          }
        

          args.log << "\nFor complete options help for a particular method, \n"
                   "use 'casm enum --desc MethodName'.\n\n";

          args.log << "Custom enumerator plugins can be added by placing source code \n"
                   "in the CASM project directory: \n"
                   "  " << primclex->dir().enumerators() << " \n\n"

                   "For examples of how to write enumerators see: \n"
                   "  $REPO/include/casm/enumerators \n"
                   "  $REPO/src/casm/enumerators \n"
                   "where: \n"
                   "  REPO=https://github.com/prisms-center/CASMcode/tree/master \n\n";
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
    if(input.size() != 1) {
      args.err_log.error("Expected a single attribute in the input");
      args.err_log << std::endl;
    }
    auto it = input.begin();
    int res = enumerators->find(it.name())->run(*primclex, *it);

    args.log << std::endl;

    return res;
  };

}

