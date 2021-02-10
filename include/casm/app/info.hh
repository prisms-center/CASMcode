#ifndef CASM_info
#define CASM_info

#include "casm/app/APICommand.hh"
#include "casm/app/info/InfoInterface.hh"
#include "casm/app/info/standard_info_method_interfaces.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

/// 'casm info' implementation
///
/// Output of information about CASM objects -- meant to be a quick way of
/// implementing new functions in C++ for use via CLI or Python, etc.
///
/// There are no enforced limits on what the methods may do, but the current
/// intention is that the input is simple JSON, with optional PrimClex pointer
/// if provided by the API call.
///
/// The input should include a list of "properties" to query via
/// DataFormatters, with default empty array indicating to print all properties
/// in a DataFormatterDictionary. For example:
/// \code
/// {
///   "prim": <prim.json format>, // optional, else current project prim
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
///
class InfoCommand : public APICommand<Completer::InfoOption> {
 public:
  static const std::string name;

  InfoCommand(const CommandArgs &_args, Completer::InfoOption &_opt);

  int vm_count_check() const override;

  int help() const override;

  int desc() const override;

  int run() const override;

  // -- custom --

  InfoInterfaceVector const &info_methods() const;

  void print_names(std::ostream &sout,
                   InfoInterfaceVector const &info_methods) const;

 private:
  mutable InfoInterfaceVector m_standard_info_methods;
  mutable InfoInterfaceVector const *m_info_method_vector;
};

}  // namespace CASM

#endif
