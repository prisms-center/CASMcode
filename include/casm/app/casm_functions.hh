#ifndef CASM_FUNCTIONS_HH
#define CASM_FUNCTIONS_HH

#include <wordexp.h>
#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/Log.hh"

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

  class PrimClex;

  /**
   *  \defgroup API
   *
   *  \brief Relates to the CASM API
   *
   *  @{
   */

  /// \brief Data structure holding basic CASM command info
  struct CommandArgs : public Logging {

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
      CommandArgs(_argc, _argv, _primclex, _root, Logging(_log, _log, _err_log)) {}


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

    int argc;
    char **argv;
    PrimClex *primclex;
    fs::path root;
    Log &log;
    Log &err_log;

    /// stores error codes when attempting to parse std::string _args -> argc, argv
    int parse_result;

    bool is_help;
    bool write_log;
    std::string command;

  private:

    void _init();

    /// Used when parsing std::string args -> argc, argv
    bool m_free_p;
    wordexp_t m_p;

  };

  typedef std::function<int (const CommandArgs &)> Command;
  typedef std::map<std::string, Command> CommandMap;

  /// \brief Return static CommandMap containing all CASM API commands
  CommandMap &command_map();

  /// \brief Executes CASM commands specified by args
  int casm_api(const CommandArgs &args);

  // /// \brief Executes casm_api in specified working directory
  // int casm_api(const CommandArgs &args, fs::path working_dir);



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

  int bset_command(const CommandArgs &args);

  int composition_command(const CommandArgs &args);

  int enum_command(const CommandArgs &args);

  int files_command(const CommandArgs &args);

  int format_command(const CommandArgs &args);

  int import_command(const CommandArgs &args);

  int init_command(const CommandArgs &args);

  int monte_command(const CommandArgs &args);

  int perturb_command(const CommandArgs &args);

  int query_command(const CommandArgs &args);

  int ref_command(const CommandArgs &args);

  int rm_command(const CommandArgs &args);

  int run_command(const CommandArgs &args);

  int select_command(const CommandArgs &args);

  int settings_command(const CommandArgs &args);

  int status_command(const CommandArgs &args);

  int super_command(const CommandArgs &args);

  int sym_command(const CommandArgs &args);

  int update_command(const CommandArgs &args);

  int version_command(const CommandArgs &args);

  int view_command(const CommandArgs &args);

  /** @} */
}

#endif
