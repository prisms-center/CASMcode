#include "casm/app/casm_functions.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigSelection.hh"
#include "casm/external/gzstream/gzstream.h"
#include "casm/version/version.hh"

namespace CASM {

  CommandArgs::CommandArgs(int _argc,
                           char *_argv[],
                           PrimClex *_primclex,
                           Log &_log,
                           Log &_err_log) :
    argc(_argc),
    argv(_argv),
    primclex(_primclex),
    log(_log),
    err_log(_err_log) {

    // set project 'root'
    root = find_casmroot(fs::current_path());

    // set 'command'
    command = (argc > 1) ? std::string(argv[1]) : "";

    // check if 'help' command
    std::vector<std::string> help_commands {
      "help",
      "-h",
      "--help",
      "version",
      "--version",
      "status",
      "format"
    };
    is_help = contains(help_commands, command);

    // check if LOG should be written
    write_log = !(is_help || root.empty() || command == "init");

  }

  /// \brief Return static CommandMap containing all CASM API commands
  CommandMap &command_map() {
    static CommandMap command_map {
      {"version", version_command},
      {"--version", version_command},
      {"status", status_command},
      {"format", format_command},
      {"init", init_command},
      {"settings", settings_command},
      {"sym", sym_command},
      {"composition", composition_command},
      {"ref", ref_command},
      {"update", update_command},
      {"enum", enum_command},
      {"super", super_command},
      {"select", select_command},
      {"bset", bset_command},
      {"perturb", perturb_command},
      {"run", run_command},
      {"query", query_command},
      {"files", files_command},
      {"import", import_command},
      {"monte", monte_command},
      {"view", view_command},
      {"help", help_command}
    };

    return command_map;
  }

  namespace api_impl {

    std::string date_time() {
      auto t = std::time(nullptr);
      char str[255];
      strftime(str, sizeof(str), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
      return std::string(str);
    };

    void write_LOG_begin(const CommandArgs &args) {

      // If not a 'version', or 'help' command, write to LOG
      fs::ofstream log(args.root / "LOG", std::ofstream::out | std::ofstream::app);
      log << "# " << date_time() << std::endl;

      // record whoami@hostname

      std::string whoami, hostname;
      {
        FILE *fp;
        char path[PATH_MAX];
        fp = popen("whoami", "r");
        while(fgets(path, PATH_MAX, fp) != NULL) {
          whoami = std::string(path);
          whoami.resize(whoami.size() - 1);
        }
      }
      {
        FILE *fp;
        char path[PATH_MAX];
        fp = popen("hostname", "r");
        while(fgets(path, PATH_MAX, fp) != NULL) {
          hostname = std::string(path);
          hostname.resize(hostname.size() - 1);
        }
      }
      log << "# " << whoami << "@" << hostname << "\n";

      // record CASM version number
      log << "# " << version()  << "\n";

      // record command arguments
      for(int i = 0; i < args.argc; i++) {
        log << args.argv[i] << " ";
      }
      log << "\n";
      log.close();
    }

    void write_LOG_end(const CommandArgs &args, int retcode) {
      fs::ofstream log(args.root / "LOG", std::ofstream::out | std::ofstream::app);
      log << "# return: " << retcode << " runtime(s): " << args.log.time_s() << "\n\n";
      log.close();
    }

  }

  /// \brief Executes CASM commands specified by args
  int casm_api(const CommandArgs &args) {

    if(args.argc == 1) {
      help_command(args);
      return ERR_INVALID_ARG;
    }

    auto it = command_map().find(args.command);
    if(it != command_map().end()) {

      if(args.write_log) {
        api_impl::write_LOG_begin(args);
      }

      args.log.restart_clock();
      int retcode = it->second(args);
      if(args.write_log) {
        api_impl::write_LOG_end(args, retcode);
      }
      return retcode;
    }
    else {
      help_command(args);
      return ERR_INVALID_ARG;
    }
  }


  /// \brief If !_primclex, construct new PrimClex stored in uniq_primclex, then
  ///        return reference to existing or constructed PrimClex
  ///
  /// \param args CommandArgs reference
  /// \param uniq_primclex Reference to null std::unique_ptr<PrimClex> to manage PrimClex, if it is constructed
  ///
  /// \returns reference to PrimClex (either newly constructed managed by
  ///          uniq_primclex, or existing pointed at by args.primclex)
  ///
  PrimClex &make_primclex_if_not(const CommandArgs &args, std::unique_ptr<PrimClex> &uniq_primclex) {
    return make_primclex_if_not(args, uniq_primclex, args.log);
  }

  /// \brief If !_primclex, construct new PrimClex stored in uniq_primclex, then
  ///        return reference to existing or constructed PrimClex
  ///
  /// \param args CommandArgs reference
  /// \param uniq_primclex Reference to null std::unique_ptr<PrimClex> to manage PrimClex, if it is constructed
  /// \param status_log where to print PrimClex construction messages
  ///
  /// \returns reference to PrimClex (either newly constructed managed by
  ///          uniq_primclex, or existing pointed at by args.primclex)
  ///
  PrimClex &make_primclex_if_not(const CommandArgs &args, std::unique_ptr<PrimClex> &uniq_primclex, Log &status_log) {
    if(!args.primclex) {
      uniq_primclex.reset(new PrimClex(args.root, status_log));
      return *uniq_primclex;
    }
    return *args.primclex;
  }

  /// \brief Return a reference to proper std::ostream
  ///
  /// \param output Output mode: False: use 'sout', True: check 'out_path' and 'gzip' to decide
  /// \param sout stream to use if not writing to file
  /// \param fout will be given an open file if writing to file
  /// \param out_path, where to write if 'output': if "STDOUT", use 'sout'; otherwise filename
  /// \param gzip: if true, write to gzip file
  ///
  /// \return reference to stream to use
  ///
  std::ostream &make_ostream_if(bool output, std::ostream &sout, std::unique_ptr<std::ostream> &fout, fs::path out_path, bool gzip) {

    if(output) {
      if(out_path.string() == "STDOUT") {
        return sout;
      }

      out_path = fs::absolute(out_path);

      if(gzip) {
        fout.reset(new gz::ogzstream(out_path.string().c_str()));
        return *fout;
      }

      fout.reset(new fs::ofstream(out_path));
      return *fout;
    }
    else {
      return sout;
    }

  }

  /// \brief Print CASM help info to args.log
  int help_command(const CommandArgs &args) {
    args.log.custom("casm usage");
    args.log << "\n";

    args.log << "casm [--version] <command> [options] [args]" << std::endl << std::endl;
    args.log << "available commands:" << std::endl;

    std::vector<std::string> subcom;
    for(auto it = command_map().begin(); it != command_map().end(); ++it) {
      std::string s = it->first;
      // check for --version
      if(s[0] != '-') {
        subcom.push_back(std::string("  ") + s);
      }
    }

    std::sort(subcom.begin(), subcom.end());
    for(auto it = subcom.begin(); it != subcom.end(); ++it) {
      args.log << *it << "\n";
    }
    args.log << "\n";

    args.log << "For help using a command: 'casm <command> --help'" << std::endl << std::endl;
    args.log << "For step by step help use: 'casm status -n'" << std::endl << std::endl;

    return 0;
  };

  int version_command(const CommandArgs &args) {
    args.log << "casm version: " << version() << std::endl;
    return 0;
  }




}
