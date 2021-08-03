#ifndef CASM_Monte2Interface
#define CASM_Monte2Interface

#include <string>

#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

class Log;
class PrimClex;
class ParentInputParser;
class jsonParser;

namespace Completer {
class Monte2Option;
}

/// Interface for `casm enum` methods
///
/// Notes:
/// - See the interfaces implemented in `casm/app/enum/methods` for examples
class Monte2InterfaceBase : public notstd::Cloneable {
  ABSTRACT_CLONEABLE(Monte2InterfaceBase)
 public:
  /// Monte Carlo method description. What will be printed by `casm enum --desc
  /// MethodName`.
  virtual std::string desc() const = 0;

  /// Monte Carlo method name (i.e. "Canonical")
  virtual std::string name() const = 0;

  /// Monte Carlo method implementation
  ///
  /// \param primclex PrimClex gives access to project resources
  /// \param json_options JSON input, as by --input or --settings CLI options
  /// \param cli_options_as_json CLI options converted to JSON. Includes:
  ///        --desc as "desc": array of string (method names to print help for)
  ///        --help as "help": bool (print/return help, including list of
  ///        available methods)
  ///        --method as "method": string (name of enumeration method)
  ///        --settings as "settings": string (path to settings JSON file)
  ///        --input as "input": string (a JSON string)
  ///        --verbosity as "verbosity": string (requested verbosity level)
  ///
  /// Notes:
  /// - It is up to the individual method to determine how to combine and use
  /// `json_options` and
  ///   `cli_options_as_json`, but by convention prefer using CLI input in case
  ///   of collisions.
  virtual void run(PrimClex &primclex, jsonParser const &json_options,
                   jsonParser const &cli_options_as_json) const = 0;
};

/// Convert `casm monte2` CLI input to JSON
jsonParser &to_json(const Completer::Monte2Option &enum_opt, jsonParser &json);

/// Combine --input / --settings JSON with CLI options
ParentInputParser make_monte2_parent_parser(
    Log &log, jsonParser const &json_options,
    jsonParser const &cli_options_as_json);

}  // namespace CASM

#endif
