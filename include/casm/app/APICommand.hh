#ifndef CASM_APICommand
#define CASM_APICommand

#include <boost/program_options.hpp>
#include "casm/app/casm_functions.hh"
#include "casm/casm_io/Log.hh"

namespace CASM {

  class PrimClex;

  /// Base class for implementing CASM API command classes
  ///
  /// - Do not inherit from this directly, instead inherit from APICommand<OptionType>
  /// - See the run_api_command function for context when implementing virtual functions
  class APICommandBase : public Logging {

  public:
    APICommandBase(const CommandArgs &_args);

    virtual ~APICommandBase() {}

    const CommandArgs &args() const;

    fs::path root() const;

    bool in_project() const;

    PrimClex &primclex() const;

    PrimClex &primclex(Log &status_log) const;

    virtual int vm_count_check() const = 0;

    virtual int help() const = 0;

    virtual int desc() const = 0;

    virtual int run() const = 0;

  private:
    const CommandArgs &m_args;
    bool m_in_project;
    mutable std::unique_ptr<PrimClex> m_primclex;
  };


  /// Base class for implementing CASM API command classes
  ///
  /// - Includes access to Option class and variables_map
  template<typename _OptionType>
  class APICommand : public APICommandBase {

  public:

    typedef _OptionType OptionType;

    APICommand(const CommandArgs &_args, OptionType &_opt) :
      APICommandBase(_args),
      m_opt(_opt) {}

    int count(std::string s) const {
      return opt().vm().count(s);
    }

    const po::variables_map &vm() const {
      return opt().vm();
    }

    const OptionType &opt() const {
      return m_opt;
    }

  private:

    OptionType &m_opt;
  };

  /// Standardizes how 'casm X' api commands are executed and implemented
  template<typename CommandType>
  int run_api_command(const CommandArgs &args) {

    typename CommandType::OptionType opt;

    try {
      po::store(po::parse_command_line(args.argc, args.argv, opt.desc()), opt.vm());

      // gets default values
      po::notify(opt.vm());

      CommandType f(args, opt);

      // checks that can be made without getting defaults
      if(!f.vm().count("help")) {
        int res = f.vm_count_check();
        if(res) {
          f.help();
          return res;
        }
      }
      // help
      else {
        return f.help();
      }

      // extended command descriptions
      if(f.vm().count("desc")) {
        return f.desc();
      }

      // main command logic:
      return f.run();

    }
    catch(po::error &e) {
      args.err_log << opt.desc() << std::endl;
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(CASM::runtime_error &e) {
      args.err_log << opt.desc() << std::endl;
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      return e.code();
    }
    catch(std::exception &e) {
      args.err_log << opt.desc() << std::endl;
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_UNKNOWN;
    }
  }

}

#endif
