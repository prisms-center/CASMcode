#ifndef CASM_FUNCTIONS_HH
#define CASM_FUNCTIONS_HH

#include <boost/filesystem/path.hpp>
#include "casm/global/definitions.hh"
#include "casm/casm_io/Log.hh"
#include "casm/app/CLIParse.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/misc/unique_cloneable_map.hh"

/**
 *  \defgroup API
 *
 *  \brief Relates to the CASM API
 *
 *  All primary CASM functionality is included in the CASM library 'libcasm'.
 *
 *  The CASM API refers to the actions that can be performed by passing string
 *  commands through the extern "C" function ::casm_capi which in turn calls
 *  functions in libcasm. The ::casm_capi function and a limited set of others
 *  that allow for constructing PrimClex objects and input/output streams are
 *  available in the library 'libccasm'.
 *
 *  The CASM API is primarily intended to be used via the command line
 *  executable 'casm', which provides documentation of the allowed options, or
 *  the 'casm' Python package, but can be also be accessed directly to enable
 *  integration with other software.
 *
 *  @{
 */

// Command line input is not valid
#define ERR_INVALID_ARG 1

// Misc. or Unknown error
#define ERR_UNKNOWN 2

// No CASM project can be found in expected location
#define ERR_NO_PROJ 3

// An expected input file is invalid
#define ERR_INVALID_INPUT_FILE 4

// An expected input file can not be found
#define ERR_MISSING_INPUT_FILE 5

// A file might be overwritten
#define ERR_EXISTING_FILE 6

// Requested command can not be performed because some dependency needs to be
//   done first (i.e. no basis set, so can't use clexulator)
#define ERR_MISSING_DEPENDS 7

// Unknown attempting to overwrite another CASM project
#define ERR_OTHER_PROJ 8

/** @} */

/// \brief Main CASM namespace
namespace CASM {

  void print_splash(std::ostream &out);

  class PrimClex;
  class jsonParser;

  class runtime_error : public std::runtime_error {
  public:
    runtime_error(std::string _what, int _code = ERR_UNKNOWN) :
      std::runtime_error(_what),
      m_code(_code) {}

    virtual ~runtime_error() {}

    int code() const {
      return m_code;
    }

  private:
    std::string m_what;
    int m_code;
  };

  /**
   *  \defgroup API
   *
   *  \brief Relates to the CASM API
   *
   *  @{
   */

  /// \brief Data structure holding basic CASM command info
  struct CommandArgs : public CLIParse {

    /// \brief CommandArgs constructor
    CommandArgs(int _argc,
                char *_argv[],
                PrimClex *_primclex,
                fs::path _root,
                const Logging &logging);

    /// \brief CommandArgs constructor
    CommandArgs(int _argc,
                char *_argv[],
                PrimClex *_primclex = nullptr,
                fs::path _root = fs::path(),
                Log &_log = default_log(),
                Log &_err_log = default_err_log()) :
      CommandArgs(_argc, _argv, _primclex, _root, Logging(_log, _log, _err_log)) {
    }


    /// \brief CommandArgs constructor
    CommandArgs(std::string _args,
                PrimClex *_primclex,
                fs::path _root,
                const Logging &logging);

    /// \brief CommandArgs constructor
    CommandArgs(std::string _args,
                PrimClex *_primclex = nullptr,
                fs::path _root = fs::path(),
                Log &_log = default_log(),
                Log &_err_log = default_err_log()) :
      CommandArgs(_args, _primclex, _root, Logging(_log, _log, _err_log)) {}

    CommandArgs(const CommandArgs &other) = delete;
    CommandArgs(CommandArgs &&other) = delete;
    CommandArgs &operator=(const CommandArgs &) = delete;
    CommandArgs &operator=(CommandArgs &&) = delete;

    /// \brief CommandArgs destructor
    ~CommandArgs();


    PrimClex *primclex;

    fs::path root;

    bool is_help;

    // if LOG should be written
    bool write_log;

    // which casm sub-command ('init', 'enum', etc.)
    std::string command;

  private:

    void _init();

  };

  typedef std::function<int (const CommandArgs &)> Command;
  typedef std::map<std::string, Command> CommandMap;

  /// \brief Return static CommandMap containing all CASM API commands
  CommandMap &command_map();

  /// \brief Executes CASM commands specified by args
  int casm_api(const CommandArgs &args);

  /// \brief If !_primclex, construct new PrimClex stored in uniq_primclex, then
  ///        return reference to existing or constructed PrimClex
  PrimClex &make_primclex_if_not(const CommandArgs &args, std::unique_ptr<PrimClex> &uniq_primclex);

  /// \brief If !_primclex, construct new PrimClex stored in uniq_primclex, then
  ///        return reference to existing or constructed PrimClex
  PrimClex &make_primclex_if_not(const CommandArgs &args, std::unique_ptr<PrimClex> &uniq_primclex, Log &status_log);

  /// \brief Return a reference to proper std::ostream
  std::ostream &make_ostream_if(
    bool output,
    std::ostream &sout,
    std::unique_ptr<std::ostream> &fout,
    fs::path out_path,
    bool gzip);


  int help_command(const CommandArgs &args);

  int version_command(const CommandArgs &args);

  template<typename OptionType>
  class InterfaceBase;

  /// \brief Used to hold a list of all algorithms that may be accessed via the API
  template<typename OptionType>
  using InterfaceMap = notstd::unique_cloneable_map<std::string, InterfaceBase<OptionType> >;

  /// \brief Base class for generic use of algorithms through the API
  template<typename OptionType>
  class InterfaceBase {

  public:

    InterfaceBase() {}

    virtual ~InterfaceBase() {}

    virtual std::string help() const = 0;

    virtual std::string name() const = 0;

    virtual int run(const PrimClex &primclex, const jsonParser &kwargs, const OptionType &opt, InterfaceMap<OptionType> const *interface_map = nullptr) const = 0;

    std::unique_ptr<InterfaceBase> clone() const {
      return std::unique_ptr<InterfaceBase>(this->_clone());
    }

  private:

    virtual InterfaceBase *_clone() const = 0;

  };

  /// \brief Use to construct an InterfaceMap
  template<typename OptionType>
  std::unique_ptr<InterfaceMap<OptionType> > make_interface_map() {

    return notstd::make_unique<InterfaceMap<OptionType> >(
    [](const InterfaceBase<OptionType> &e) -> std::string {
      return e.name();
    },
    [](const InterfaceBase<OptionType> &e) -> notstd::cloneable_ptr<InterfaceBase<OptionType> > {
      return notstd::clone(e);
    });
  }

  /** @} */
}

#endif
