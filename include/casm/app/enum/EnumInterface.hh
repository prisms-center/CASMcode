#ifndef CASM_EnumInterface
#define CASM_EnumInterface

#include <string>
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

  class Log;
  class PrimClex;
  class ParentInputParser;
  class jsonParser;

  namespace Completer {
    class EnumOption;
  }

  /// Interface for `casm enum` methods
  ///
  /// Notes:
  /// - See the interfaces implemented in `casm/app/enum/methods` for examples
  class EnumInterfaceBase : public notstd::Cloneable {
    ABSTRACT_CLONEABLE(EnumInterfaceBase)
  public:

    /// Enumeration method description. What will be printed by `casm enum --desc MethodName`.
    virtual std::string desc() const = 0;

    /// Enumeration method name (i.e. "ConfigEnumAllOccupations")
    virtual std::string name() const = 0;

    /// Enumeration method implementation
    ///
    /// \param primclex PrimClex gives access to project resources
    /// \param json_options JSON input, as by --input or --settings CLI options
    /// \param cli_options_as_json CLI options converted to JSON. Includes:
    ///        --desc as "desc": array of string (method names to print help for)
    ///        --help as "help": bool (print/return help, including list of available methods)
    ///        --method as "method": string (name of enumeration method)
    ///        --settings as "settings": string (path to settings JSON file)
    ///        --input as "input": string (a JSON string)
    ///        --min as "min": int (meaning depends on method, if any)
    ///        --max as "max": int (meaning depends on method, if any)
    ///        --all as "all": bool (meaning depends on method, if any)
    ///        --scelnames as "scelnames": array of string (meaning depends on method, if any)
    ///        --confignames as "confignames": array of string (meaning depends on method, if any)
    ///        --verbosity as "verbosity": string (requested verbosity level)
    ///        --filter as "filter": array of string (query expression used filter results)
    ///        --dry-run as "dry_run": bool (if true, do not commit database upon completion)
    ///
    /// Notes:
    /// - It is up to the individual method to determine how to combine and use `json_options` and
    ///   `cli_options_as_json`, but by convention prefer using CLI input in case of collisions.
    virtual void run(PrimClex &primclex,
                     jsonParser const &json_options,
                     jsonParser const &cli_options_as_json) const = 0;

  };

  /// Convert `casm enum` CLI input to JSON
  jsonParser &to_json(const Completer::EnumOption &enum_opt, jsonParser &json);

  /// Combine --input / --settings JSON with CLI options
  ParentInputParser make_enum_parent_parser(Log &log,
                                            jsonParser const &json_options,
                                            jsonParser const &cli_options_as_json);

  /// Standardized string that can prepended to method output to make clear that method results are
  /// not going to be committed
  inline std::string dry_run_msg(bool dry_run) {
    return dry_run ? "(dry run) " : "";
  }

}

#endif
