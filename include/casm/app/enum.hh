#ifndef CASM_enum
#define CASM_enum

#include "casm/app/APICommand.hh"
#include "casm/app/enum/EnumInterface.hh"
#include "casm/app/enum/standard_enumerator_interfaces.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

  /// 'casm enum' implementation
  ///
  /// The `enum` command provides access to methods that construct and save Supercells or
  /// Configurations.
  ///
  /// New "standard" methods may be implemented by:
  /// - Implement an interface derived from `EnumInteraceBase`. This involves:
  ///   - Giving the method a "name"
  ///   - Provide a "help" string describing the method and JSON input format and options
  ///   - Implement a "run" method that accepts inputs and executes the method
  ///   - See the interfaces implemented in `casm/app/enum/methods` for examples
  /// - Insert the interface into the standard enumerator map by adding it to the implementation
  ///   of `make_standard_enumerator_interfaces`.
  ///
  /// Custom enumeration methods can also be implemented via plugins. TODO: Full description.
  class EnumCommand : public APICommand<Completer::EnumOption> {

  public:

    static const std::string name;

    EnumCommand(CommandArgs const &_args, Completer::EnumOption &_opt);

    int vm_count_check() const override;

    int help() const override;

    int desc() const override;

    int run() const override;

    // -- custom --

    EnumInterfaceVector const &enumerators() const;

    void print_names(std::ostream &sout, EnumInterfaceVector const &enumerators) const;

  private:

    mutable EnumInterfaceVector m_standard_enumerators;
    mutable EnumInterfaceVector const *m_enumerator_vector;
  };
}

#endif
