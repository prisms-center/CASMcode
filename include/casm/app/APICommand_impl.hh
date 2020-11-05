#ifndef CASM_APICommand_impl
#define CASM_APICommand_impl

#include "casm/app/APICommand.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {

  /// Parse command line options and make API command. Throws for parsing errors.
  template<typename CommandType>
  std::unique_ptr<CommandType> make_api_command(CommandArgs const &args, typename CommandType::OptionType &opt) {
    po::store(po::parse_command_line(args.argc(), args.argv(), opt.desc()), opt.vm());
    po::notify(opt.vm());
    return notstd::make_unique<CommandType>(args, opt);
  }

  /// Parse command line options and make API command. Throws for parsing errors.
  ///
  /// \param cli_str CLI args string, ex: 'casm X ...'
  /// \param primclex Existing PrimClex
  ///
  template<typename CommandType>
  std::unique_ptr<CommandType> make_api_command(std::string cli_str, PrimClex &primclex, typename CommandType::OptionType &opt) {
    fs::path root;
    if(primclex.has_dir()) {
      root = primclex.dir().root_dir();
    }
    CommandArgs args {cli_str, &primclex, root};
    return make_api_command<CommandType>(args, opt);
  }

  /// Standardizes how 'casm X' api commands are executed and implemented
  template<typename CommandType>
  int run_api_command(const CommandArgs &args) {

    typename CommandType::OptionType opt;

    try {
      po::store(po::parse_command_line(args.argc(), args.argv(), opt.desc()), opt.vm());

      // gets default values
      po::notify(opt.vm());

      CommandType f {args, opt};

      int code;

      // help
      if(f.vm().count("help")) {
        code = f.help();
      }
      // extended command descriptions
      else if(f.vm().count("desc")) {
        code = f.desc();
      }
      // checks that can be made without getting defaults
      else {
        code = f.vm_count_check();
        if(code) {
          f.help();
        }
        else {
          code = f.run();
        }
      }

      return code;
    }
    catch(po::error &e) {
      err_log() << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(CASM::runtime_error &e) {
      err_log() << "ERROR: " << e.what() << std::endl << std::endl;
      return e.code();
    }
    catch(std::exception &e) {
      err_log() << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_UNKNOWN;
    }
  }

}

#endif
