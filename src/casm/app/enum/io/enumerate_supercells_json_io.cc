#include "casm/app/enum/io/enumerate_supercells_json_io.hh"

#include "casm/app/enum/enumerate_supercells.hh"
#include "casm/casm_io/dataformatter/DataFormatterFilter_impl.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/dataformatter/FormattedDataFile_impl.hh"
#include "casm/casm_io/json/InputParser_impl.hh"

namespace CASM {

/// A standard approach to combine CLI options with user JSON options into one
/// JSON document
///
/// Not all supercell enumerator interfaces need to incorporate the
/// CLI output using this method, but it is the most common way to do so (only
/// way currently).
jsonParser combine_supercell_enum_json_options(
    jsonParser const &json_options, jsonParser const &cli_options_as_json) {
  // combine JSON options and CLI options
  jsonParser json_combined{json_options};

  if (cli_options_as_json.contains("min")) {
    json_combined["min"] = cli_options_as_json["min"];
  }
  if (cli_options_as_json.contains("max")) {
    json_combined["max"] = cli_options_as_json["max"];
  }
  if (cli_options_as_json.contains("filter")) {
    json_combined["filter"] = cli_options_as_json["filter"];
  }
  if (cli_options_as_json.contains("dry_run")) {
    json_combined["dry_run"] = cli_options_as_json["dry_run"];
  }
  if (cli_options_as_json.contains("verbosity")) {
    json_combined["verbosity"] = cli_options_as_json["verbosity"];
  }

  return json_combined;
}

// Enable InputParser<EnumerateSupercellsOptions>
void parse(InputParser<EnumerateSupercellsOptions> &parser,
           std::string method_name, PrimClex const &primclex,
           DataFormatterDictionary<Supercell> const &dict) {
  parser.value = notstd::make_unique<EnumerateSupercellsOptions>(primclex);
  auto &options = *parser.value;

  options.method_name = method_name;

  parser.optional_else(options.dry_run, "dry_run", false);

  options.verbosity = parse_verbosity(parser);

  std::vector<std::string> filter_expression;
  parser.optional(filter_expression, "filter");
  if (filter_expression.size()) {
    options.filter = make_data_formatter_filter(filter_expression, dict);
  }

  parser.optional_else(options.output_supercells, "output_supercells", false);

  if (options.output_supercells) {
    // TODO: separate parser for FormattedDataFileOptions

    fs::path base{"output_supercells_options"};

    std::string file_path_str;
    parser.optional_else<std::string>(file_path_str, base / "path", "enum.out");
    fs::path file_path{file_path_str};

    bool json_output;
    parser.optional_else(json_output, base / "json", false);

    bool json_arrays_output;
    parser.optional_else(json_arrays_output, base / "json_arrays", false);

    bool compress;
    parser.optional_else(compress, base / "compress", false);

    if (compress) {
      if (file_path.extension() != ".gz" && file_path.extension() != ".GZ") {
        file_path += ".gz";
      }
    }

    bool output_filtered_supercells;
    parser.optional_else(output_filtered_supercells,
                         base / "output_filtered_supercells", false);

    options.output_options = FormattedDataFileOptions{
        file_path, json_output, json_arrays_output, compress};
    options.output_filtered_supercells = output_filtered_supercells;
  }
}

}  // namespace CASM
