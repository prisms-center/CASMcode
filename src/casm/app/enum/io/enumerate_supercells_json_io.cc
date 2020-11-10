#include "casm/app/enum/enumerate_supercells.hh"
#include "casm/app/enum/io/enumerate_supercells_json_io.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/dataformatter/DataFormatterFilter_impl.hh"
#include "casm/casm_io/json/InputParser_impl.hh"

namespace CASM {

  /// A standard approach to combine CLI options with user JSON options into one JSON document
  ///
  /// Not all supercell enumerator interfaces need to incorporate the
  /// CLI output using this method, but it is the most common way to do so (only way currently).
  jsonParser combine_supercell_enum_json_options(
    jsonParser const &json_options,
    jsonParser const &cli_options_as_json) {

    // combine JSON options and CLI options
    jsonParser json_combined {json_options};

    if(cli_options_as_json.contains("min")) {
      json_combined["min"] = cli_options_as_json["min"];
    }
    if(cli_options_as_json.contains("max")) {
      json_combined["max"] = cli_options_as_json["max"];
    }
    if(cli_options_as_json.contains("filter")) {
      json_combined["filter"] = cli_options_as_json["filter"];
    }
    if(cli_options_as_json.contains("dry_run")) {
      json_combined["dry_run"] = cli_options_as_json["dry_run"];
    }
    if(cli_options_as_json.contains("verbosity")) {
      json_combined["verbosity"] = cli_options_as_json["verbosity"];
    }

    return json_combined;
  }

  // Enable InputParser<EnumerateSupercellsOptions>
  void parse(
    InputParser<EnumerateSupercellsOptions> &parser,
    std::string method_name,
    PrimClex const &primclex,
    DataFormatterDictionary<Supercell> const &dict) {

    parser.value = notstd::make_unique<EnumerateSupercellsOptions>(primclex);
    auto &options = *parser.value;

    options.method_name = method_name;

    parser.optional_else(options.dry_run, "dry_run", false);

    options.verbosity = parse_verbosity(parser);

    std::vector<std::string> filter_expression;
    parser.optional(filter_expression, "filter");
    if(filter_expression.size()) {
      options.filter = make_data_formatter_filter(filter_expression, dict);
    }
  }

}
