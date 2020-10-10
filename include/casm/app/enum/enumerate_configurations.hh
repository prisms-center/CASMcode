#ifndef CASM_enum_enumerate_configurations
#define CASM_enum_enumerate_configurations

#include <functional>
#include <map>
#include <string>

namespace CASM {

  class Configuration;
  class PrimClex;
  class Supercell;

  namespace DB {
    template<typename T> class Database;
  }

  /// Options for the `enumerate_configurations` function
  struct EnumerateConfigurationsOptions {

    EnumerateConfigurationsOptions(PrimClex const &primclex):
      primclex_ptr(&primclex) {}

    /// Method name, for printing progress
    std::string method_name;

    /// If `primitive_only==true`, only the primitive configuration is inserted,
    /// otherwise non-primitive are also inserted
    bool primitive_only = true;

    /// If `filter(configuration)==true`, keep configuration, else skip
    std::function<bool (Configuration const &)> filter;

    /// If `dry_run==true`, do not save results, just print to screen
    bool dry_run = false;

    /// Printing verbosity level
    int verbosity = 10;

    /// Use while transitioning Supercell to no longer need a `PrimClex const *`
    PrimClex const *primclex_ptr = nullptr;

  };

  /// Enumerate configurations
  template<typename MakeEnumeratorFunctor, typename InputType>
  void enumerate_configurations(
    EnumerateConfigurationsOptions const &options,
    MakeEnumeratorFunctor make_enumerator_f,
    std::map<std::string, InputType> input_name_value_map,
    DB::Database<Supercell> &supercell_db,
    DB::Database<Configuration> &configuration_db);

}

#endif
