#ifndef CASM_FUNCTIONS_HH
#define CASM_FUNCTIONS_HH

#include <boost/filesystem/path.hpp>

#include "casm/app/CLIParse.hh"
#include "casm/app/errors.hh"
#include "casm/global/definitions.hh"
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

/** @} */

/// \brief Main CASM namespace
namespace CASM {

void print_splash(std::ostream &out);

class PrimClex;
class jsonParser;

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
  CommandArgs(int _argc, char *_argv[], PrimClex *_primclex = nullptr,
              fs::path _root = fs::path());

  /// \brief CommandArgs constructor
  CommandArgs(std::string _args, PrimClex *_primclex = nullptr,
              fs::path _root = fs::path());

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

typedef std::function<int(const CommandArgs &)> Command;
typedef std::map<std::string, Command> CommandMap;

/// \brief Return static CommandMap containing all CASM API commands
CommandMap &command_map();

/// \brief Executes CASM commands specified by args
int casm_api(const CommandArgs &args);

/// \brief If !_primclex, construct new PrimClex stored in uniq_primclex, then
///        return reference to existing or constructed PrimClex
PrimClex &make_primclex_if_not(const CommandArgs &args,
                               std::unique_ptr<PrimClex> &uniq_primclex);

/// \brief Return a reference to proper std::ostream
std::ostream &make_ostream_if(bool output, std::ostream &sout,
                              std::unique_ptr<std::ostream> &fout,
                              fs::path out_path, bool gzip);

int help_command(const CommandArgs &args);

int version_command(const CommandArgs &args);

template <typename OptionType>
class InterfaceBase;

/// \brief Used to hold a list of all algorithms that may be accessed via the
/// API
template <typename OptionType>
using InterfaceMap =
    notstd::unique_cloneable_map<std::string, InterfaceBase<OptionType> >;

/// \brief Base class for generic use of algorithms through the API
template <typename OptionType>
class InterfaceBase {
 public:
  InterfaceBase() {}

  virtual ~InterfaceBase() {}

  virtual std::string help() const = 0;

  virtual std::string name() const = 0;

  virtual int run(PrimClex const &primclex,
                  jsonParser const &json_options) const = 0;

  std::unique_ptr<InterfaceBase> clone() const {
    return std::unique_ptr<InterfaceBase>(this->_clone());
  }

 private:
  virtual InterfaceBase *_clone() const = 0;
};

/// \brief Use to construct an InterfaceMap
template <typename OptionType>
std::unique_ptr<InterfaceMap<OptionType> > make_interface_map() {
  return notstd::make_unique<InterfaceMap<OptionType> >(
      [](const InterfaceBase<OptionType> &e) -> std::string {
        return e.name();
      },
      [](const InterfaceBase<OptionType> &e)
          -> notstd::cloneable_ptr<InterfaceBase<OptionType> > {
        return notstd::clone(e);
      });
}

/** @} */
}  // namespace CASM

#endif
