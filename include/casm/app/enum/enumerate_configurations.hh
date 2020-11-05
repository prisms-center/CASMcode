#ifndef CASM_enum_enumerate_configurations
#define CASM_enum_enumerate_configurations

#include <functional>
#include <map>
#include <string>

#include "casm/casm_io/dataformatter/FormattedDataFile.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {

  class Configuration;
  class PrimClex;
  class Supercell;

  namespace DB {
    template<typename T> class Database;
  }

  /// Options for the `enumerate_configurations` function
  struct ConfigEnumOptions {

    ConfigEnumOptions(PrimClex const &primclex):
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

    /// If true, output a selection file with information about enumerated configurations
    bool output_configurations = false;

    /// Options for construcing FormattedDataFile object
    FormattedDataFileOptions output_options;

    /// If true, include output for configurations that were filtered out
    bool output_filtered_configurations = false;

  };

  /// Collect information during `enumerate_configurations` function for optional output
  template<typename EnumeratorType, typename InitialStateType>
  struct ConfigEnumData {

    ConfigEnumData(PrimClex const &_primclex,
                   Index _initial_state_index,
                   std::string const &_initial_state_name,
                   InitialStateType const &_initial_state,
                   EnumeratorType const &_enumerator,
                   Configuration const &_configuration):
      primclex(_primclex),
      initial_state_index(_initial_state_index),
      initial_state_name(_initial_state_name),
      initial_state(_initial_state),
      enumerator(_enumerator),
      configuration(_configuration) {}

    PrimClex const &primclex;
    Index initial_state_index;
    std::string const &initial_state_name;
    InitialStateType const &initial_state;
    EnumeratorType const &enumerator;
    Configuration const &configuration;

    bool is_excluded_by_filter = false;
    ConfigInsertResult insert_result;
  };

  /// Enumerate configurations
  template <
    typename MakeEnumeratorFunction,
    typename InputNameValuePairIterator,
    typename ConfigEnumDataType >
  void enumerate_configurations(
    PrimClex const &primclex,
    ConfigEnumOptions const &options,
    MakeEnumeratorFunction make_enumerator_f,
    InputNameValuePairIterator name_value_pairs_begin,
    InputNameValuePairIterator name_value_pairs_end,
    DataFormatter<ConfigEnumDataType> const &formatter);

}

#endif
