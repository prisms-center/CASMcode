#ifndef CASM_InfoInterface
#define CASM_InfoInterface

#include <string>

#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

class PrimClex;
class jsonParser;

/// Interface for `casm enum` methods
///
/// Notes:
/// - See the interfaces implemented in `casm/app/enum/methods` for examples
class InfoInterfaceBase : public notstd::Cloneable {
  ABSTRACT_CLONEABLE(InfoInterfaceBase)
 public:
  /// Enumeration method description. What will be printed by `casm info --desc
  /// MethodName`.
  virtual std::string desc() const = 0;

  /// Enumeration method name (i.e. "prim", "supercell", "dof_space", etc.)
  virtual std::string name() const = 0;

  /// Enumeration method implementation
  ///
  /// \param json_options JSON input, as by --input or --settings CLI options
  /// \param primclex Optional PrimClex pointer gives access to project
  ///     resources. Will only by provided if provided to API, for use of
  ///     current working directory project, must be constructed by the method.
  ///
  /// It is up to the individual method to determine how to use `json_options`
  /// and document itself via `desc()`. For most methods, the input should
  /// include a list of "properties" to query via DataFormatters, with default
  /// empty array indicating to print all properties in a
  /// DataFormatterDictionary. For example:
  /// \code
  /// {
  ///   "prim": <prim.json format>, // optional
  ///   "something_else": ..., // for example Configuration JSON or other inputs
  ///   // optional, default=[] prints all DataFormatter properties///
  ///   "properties": ["name", "name", "name(args, ...)"]
  /// }
  /// \endcode
  ///
  /// The output is primarily intended to be JSON directly to the `CASM::log()`,
  /// but it is up to the each method to document itself.
  ///
  /// Implement a new `casm info` method by following the example of
  /// `PrimInfoInterface`.
  virtual void run(jsonParser const &json_options,
                   PrimClex const *primclex = nullptr) const = 0;
};

}  // namespace CASM

#endif
