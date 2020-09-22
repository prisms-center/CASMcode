#ifndef CASM_app_enum_enumerate_configurations_json_io
#define CASM_app_enum_enumerate_configurations_json_io

#include <string>
#include "casm/casm_io/dataformatter/DataFormatterDecl.hh"

namespace CASM {

  class Configuration;
  struct EnumerateConfigurationsOptions;
  template<typename T> class InputParser;
  class PrimClex;
  class jsonParser;

  /// A standard approach to combine CLI options with user JSON options into one JSON document
  jsonParser combine_configuration_enum_json_options(
    jsonParser const &json_options,
    jsonParser const &cli_options_as_json);

  // Enable InputParser<EnumerateConfigurationsOptions>
  void parse(
    InputParser<EnumerateConfigurationsOptions> &parser,
    std::string method_name,
    PrimClex const &primclex,
    DataFormatterDictionary<Configuration> const &dict);

}

#endif
