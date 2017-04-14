#include<cstring>
#include<unistd.h>

#include "casm/app/casm_functions.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/ConfigSelection.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

  namespace Completer {

    RunOption::RunOption(): OptionHandlerBase("run") {}

    void RunOption::initialize() {
      add_help_suboption();
      add_configlist_suboption();

      m_desc.add_options()
      ("exec,e", po::value<std::string>(&m_exec_str)->required()->value_name(ArgHandler::command()), "Command to execute");
      return;
    }

    const std::string &RunOption::exec_str() const {
      return m_exec_str;
    };

  }

  // ///////////////////////////////////////
  // 'run' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int run_command(const CommandArgs &args) {
    std::string exec;
    fs::path selection;
    po::variables_map vm;

    /// Set command line options using boost program_options
    Completer::RunOption run_opt;

    try {
      po::store(po::parse_command_line(args.argc, args.argv, run_opt.desc()), vm); // can throw

      /** --help option
       */
      if(vm.count("help")) {
        args.log << "\n";
        args.log << run_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        args.log << "\n";
        args.log << run_opt.desc() << std::endl;
        args.log << "DESCRIPTION\n"
                 << "    Executes the requested command for each selected configuration,\n"
                 << "    with the path to the configuration as an argument.             \n\n"
                 << "    Example: casm run --exec \"vasp.relax\"\n"
                 << "    - calls:\n"
                 << "        'vasp.relax $ROOT/training_data/$SCELNAME/$CONFIGID'\n"
                 << "      for each config selected in config_list\n\n";

        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems

      selection = run_opt.selection_path();
    }
    catch(po::error &e) {
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      args.err_log << run_opt.desc() << std::endl;
      return 1;
    }
    catch(std::exception &e) {
      args.err_log << "ERROR: " << e.what() << ".\n       Exiting..." << std::endl;
      return 1;

    }

    const fs::path &root = args.root;
    if(root.empty()) {
      args.err_log.error("No casm project found");
      args.err_log << std::endl;
      return ERR_NO_PROJ;
    }

    // If 'args.primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);
    const auto &dir = primclex.dir();

    try {
      if(!vm.count("config") || (selection == "MASTER")) {
        for(auto it = primclex.selected_config_begin(); it != primclex.selected_config_end(); ++it) {
          if(!fs::exists(dir.POS(it->name()))) {
            it->write_pos();
          }

          Popen process;

          process.popen(run_opt.exec_str() + " " + it->get_path().string());

          process.print(args.log);
        }
      }
      else if(vm.count("config") && fs::exists(selection)) {
        ConfigSelection<true> config_select(primclex, selection);
        for(auto it = config_select.selected_config_begin(); it != config_select.selected_config_end(); ++it) {
          if(!fs::exists(dir.POS(it->name()))) {
            it->write_pos();
          }

          Popen process;

          process.popen(run_opt.exec_str() + " " + it->get_path().string());

          process.print(args.log);
        }

      }
      else {
        args.err_log << "ERROR: Invalid input. Option '--config' accepts one argument (either 'MASTER' or a path to a valid configuration selection file)." << std::endl
                     << "       Exiting...\n";
        return 1;
      }
    }
    catch(std::exception &e) {
      args.err_log << "ERROR: Invalid input. Option '--config' accepts one argument (either 'MASTER' or a path to a valid configuration selection file)." << std::endl
                   << "       Exiting...\n";
      return 1;
    }




    args.log << "\n***************************\n" << std::endl;

    args.log << std::endl;

    return 0;
  };

}


