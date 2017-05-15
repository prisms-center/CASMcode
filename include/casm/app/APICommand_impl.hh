#ifndef CASM_APICommand_impl
#define CASM_APICommand_impl

#include "casm/app/APICommand.hh"

namespace CASM {

  /// Standardizes how 'casm X' api commands are executed and implemented
  template<typename CommandType>
  int run_api_command(const CommandArgs &args) {

    typename CommandType::OptionType opt;

    try {
      po::store(po::parse_command_line(args.argc(), args.argv(), opt.desc()), opt.vm());

      // gets default values
      po::notify(opt.vm());

      CommandType f(args, opt);

      // checks that can be made without getting defaults
      if(!f.vm().count("help") && !f.vm().count("desc")) {
        int res = f.vm_count_check();
        if(res) {
          f.help();
          return res;
        }
      }
      // help
      else if(f.vm().count("help")) {
        return f.help();
      }
      // extended command descriptions
      else if(f.vm().count("desc")) {
        return f.desc();
      }

      // main command logic:
      return f.run();

    }
    catch(po::error &e) {
      args.err_log() << opt.desc() << std::endl;
      args.err_log() << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(CASM::runtime_error &e) {
      args.err_log() << opt.desc() << std::endl;
      args.err_log() << "ERROR: " << e.what() << std::endl << std::endl;
      return e.code();
    }
    catch(std::exception &e) {
      args.err_log() << opt.desc() << std::endl;
      args.err_log() << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_UNKNOWN;
    }
  }

}

#endif
