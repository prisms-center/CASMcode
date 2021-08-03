#ifndef CASM_monte2
#define CASM_monte2

#include "casm/app/APICommand.hh"
#include "casm/app/monte2/Monte2Interface.hh"
#include "casm/app/monte2/standard_monte2_method_interfaces.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

/// 'casm monte2' implementation
///
/// The `monte2` command provides access to Monte Carlo methods
///
/// New "standard" methods may be implemented by:
/// - Implement an interface derived from `Monte2InteraceBase`. This involves:
///   - Giving the method a "name"
///   - Provide a "help" string describing the method and JSON input format and
///   options
///   - Implement a "run" method that accepts inputs and executes the method
///   - See the interfaces implemented in `casm/app/monte2/methods` for examples
/// - Insert the interface into the standard method map by adding it to the
/// implementation
///   of `make_standard_monte2_interfaces`.
///
/// Custom Monte Carlo methods can also be implemented via plugins. TODO: Full
/// description.
class Monte2Command : public APICommand<Completer::Monte2Option> {
 public:
  static const std::string name;

  Monte2Command(CommandArgs const &_args, Completer::Monte2Option &_opt);

  int vm_count_check() const override;

  int help() const override;

  int desc() const override;

  int run() const override;

  // -- custom --

  Monte2InterfaceVector const &methods() const;

  void print_names(std::ostream &sout,
                   Monte2InterfaceVector const &enumerators) const;

 private:
  mutable Monte2InterfaceVector m_standard_methods;
  mutable Monte2InterfaceVector const *m_method_vector;
};
}  // namespace CASM

#endif
