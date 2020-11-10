#include "casm/app/enum/enumerate_configurations.hh"
#include "casm/app/enum/io/enumerate_configurations_json_io.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/dataformatter/DataFormatterFilter_impl.hh"
#include "casm/casm_io/dataformatter/FormattedDataFile_impl.hh"
#include "casm/casm_io/json/InputParser_impl.hh"

namespace CASM {

  /// A standard approach to combine CLI options with user JSON options into one JSON document
  ///
  /// Not all configuration enumerator interfaces need to incorporate the
  /// CLI output using this method, but it is the most common way to do so.
  jsonParser combine_configuration_enum_json_options(
    jsonParser const &json_options,
    jsonParser const &cli_options_as_json) {

    jsonParser json_combined {json_options};

    ///   "min": <int, min supercell volume to do enumerations>,
    ///   "max": <int, max supercell volume to do enumerations>,
    ///   "scelnames": <array of string, list of supercell names, context dependent usage>,
    ///   "confignames": <array of string, list of config names, context dependent usage>,
    ///   "all": <bool, use all existing supercells for enumeration, ignore other inputs>,
    ///   "filter": <array of string, filter query/select expression, save configs if evaluates true>,
    ///   "verbosity": <string, to be read by Log::verbosity_level>,
    ///   "dry_run": <bool, print/return method results but do not save results>

    if(cli_options_as_json.contains("min")) {
      json_combined["supercells"]["min"] = cli_options_as_json["min"];
    }
    if(cli_options_as_json.contains("max")) {
      json_combined["supercells"]["max"] = cli_options_as_json["max"];
    }
    if(cli_options_as_json.contains("scelnames")) {
      json_combined["scelnames"] = cli_options_as_json["scelnames"];
    }
    if(cli_options_as_json.contains("confignames")) {
      json_combined["confignames"] = cli_options_as_json["scelnames"];
    }
    if(cli_options_as_json.contains("all") && cli_options_as_json["all"] == true) {
      if(json_combined.contains("supercells")) {
        json_combined.erase("supercells");
      }
      if(json_combined.contains("scelnames")) {
        json_combined.erase("scelnames");
      }
      if(json_combined.contains("confignames")) {
        json_combined.erase("confignames");
      }
      if(json_combined.contains("config_selection")) {
        json_combined.erase("config_selection");
      }
      json_combined["supercell_selection"] = "ALL";
    }
    if(cli_options_as_json.contains("filter")) {
      json_combined["filter"] = cli_options_as_json["filter"];
    }
    if(cli_options_as_json.contains("verbosity")) {
      json_combined["verbosity"] = cli_options_as_json["verbosity"];
    }
    if(cli_options_as_json.contains("dry_run")) {
      json_combined["dry_run"] = cli_options_as_json["dry_run"];
    }

    return json_combined;
  }

  // Enable InputParser<ConfigEnumOptions>
  void parse(
    InputParser<ConfigEnumOptions> &parser,
    std::string method_name,
    PrimClex const &primclex,
    DataFormatterDictionary<Configuration> const &dict) {

    parser.value = notstd::make_unique<ConfigEnumOptions>(primclex);
    auto &options = *parser.value;

    options.method_name = method_name;

    parser.optional_else(options.primitive_only, "primitive_only", true);

    parser.optional_else(options.dry_run, "dry_run", false);

    options.verbosity = parse_verbosity(parser);

    parser.optional(options.filter_expression, "filter");
    if(!options.filter_expression.empty()) {
      options.filter = make_data_formatter_filter(options.filter_expression, dict);
    }

    parser.optional_else(options.output_configurations, "output_configurations", false);

    if(options.output_configurations) {

      // TODO: separate parser for FormattedDataFileOptions

      fs::path base {"output_configurations_options"};

      std::string file_path_str;
      parser.optional_else<std::string>(file_path_str, base / "path", "enum.out");
      fs::path file_path {file_path_str};

      bool json_output;
      parser.optional_else(json_output, base / "json", false);

      bool json_arrays_output;
      parser.optional_else(json_arrays_output, base / "json_arrays", false);

      bool compress;
      parser.optional_else(compress, base / "compress", false);

      if(compress) {
        if(file_path.extension() != ".gz" && file_path.extension() != ".GZ") {
          file_path += ".gz";
        }
      }

      bool output_filtered_configurations;
      parser.optional_else(output_filtered_configurations, base / "output_filtered_configurations", false);

      options.output_options = FormattedDataFileOptions {file_path, json_output, json_arrays_output, compress};
      options.output_filtered_configurations = output_filtered_configurations;
    }
  }
}
