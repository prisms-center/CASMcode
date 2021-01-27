#ifndef CASM_app_enum_enumerate_supercells_json_io
#define CASM_app_enum_enumerate_supercells_json_io

#include <string>

#include "casm/casm_io/dataformatter/DataFormatterDecl.hh"

namespace CASM {

template <typename T>
class InputParser;
struct EnumerateSupercellsOptions;
class jsonParser;

/// A standard approach to combine CLI options with user JSON options into one
/// JSON document
jsonParser combine_supercell_enum_json_options(
    jsonParser const &json_options, jsonParser const &cli_options_as_json);

// Enable InputParser<EnumerateSupercellsOptions>
void parse(InputParser<EnumerateSupercellsOptions> &parser,
           std::string method_name, PrimClex const &primclex,
           DataFormatterDictionary<Supercell> const &dict);

}  // namespace CASM

#endif
