#include "casm/app/enum/enumerate_configurations_impl.hh"
#include "casm/app/enum/methods/ConfigEnumRandomOccupationsInterface.hh"
#include "casm/clex/ConfigEnumRandomOccupations.hh"

#include "casm/app/APICommand.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler_impl.hh"
#include "casm/app/enum/standard_ConfigEnumInput_help.hh"
#include "casm/app/enum/io/enumerate_configurations_json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"

namespace CASM {

  std::string ConfigEnumRandomOccupationsInterface::desc() const {

    std::string custom_options =
      "  n_config: integer (optional, default=100) \n"
      "    How many random configurations to generate. Includes duplicate and pre-\n"
      "    existing configurations.                                                 \n\n";

    std::string examples =
      "  Examples:\n"
      "    To enumerate 200 random occupations in supercells up to and including size 4:\n"
      "      casm enum --method ConfigEnumRandomOccupations -i \n"
      "        '{\"supercells\":{\"max\":4}, \"n_config\": 200}' \n"
      "\n"
      "    To enumerate 200 random occupations in all existing supercells:\n"
      "      casm enum --method ConfigEnumRandomOccupations --all -i '{\"n_config\": 200}' \n"
      "\n"
      "    To enumerate 100 random occupations in all existing supercells:\n"
      "      casm enum --method ConfigEnumRandomOccupations --all \n"
      "\n"
      "    To enumerate 200 random occupations in particular supercells:\n"
      "      casm enum --method ConfigEnumRandomOccupations -i \n"
      "      '{ \n"
      "        \"scelnames\": [\n"
      "           \"SCEL1_1_1_1_0_0_0\",\n"
      "           \"SCEL2_1_2_1_0_0_0\",\n"
      "           \"SCEL4_1_4_1_0_0_0\"\n"
      "        ],\n"
      "        \"n_config\": 200\n"
      "      }' \n\n";

    return name() + ": \n\n" + custom_options + standard_ConfigEnumInput_help() + examples;
  }

  std::string ConfigEnumRandomOccupationsInterface::name() const {
    return ConfigEnumRandomOccupations::enumerator_name;
  }

  void parse(InputParser<ConfigEnumRandomOccupationsParams> &parser, MTRand &mtrand) {
    Index n_config;
    parser.optional_else(n_config, "n_config", Index {100});

    parser.value = notstd::make_unique<ConfigEnumRandomOccupationsParams>(mtrand, n_config);

  }

  void ConfigEnumRandomOccupationsInterface::run(
    APICommandBase const &cmd,
    jsonParser const &json_options,
    jsonParser const &cli_options_as_json) const {

    // Get project and log
    PrimClex &primclex = cmd.primclex();
    Logging const &logging = cmd;

    // combine JSON options and CLI options
    jsonParser json_combined = combine_configuration_enum_json_options(
                                 json_options,
                                 cli_options_as_json);

    // Read input data from JSON
    ParentInputParser parser {json_combined};

    MTRand mtrand;
    auto random_occ_params_parser_ptr = parser.parse_as<ConfigEnumRandomOccupationsParams>(mtrand);
    auto input_parser_ptr = parser.parse_as<std::vector<std::pair<std::string, ConfigEnumInput>>>(
                              primclex.shared_prim(),
                              &primclex,
                              primclex.db<Supercell>(),
                              primclex.db<Configuration>());
    auto options_parser_ptr = parser.parse_as<EnumerateConfigurationsOptions>(
                                ConfigEnumRandomOccupations::enumerator_name,
                                primclex,
                                primclex.settings().query_handler<Configuration>().dict());

    std::runtime_error error_if_invalid {"Error reading ConfigEnumRandomOccupations JSON input"};
    report_and_throw_if_invalid(parser, logging.log(), error_if_invalid);

    ConfigEnumRandomOccupationsParams const &random_occ_params = *random_occ_params_parser_ptr->value;
    auto const &input_name_value_pairs = *input_parser_ptr->value;
    EnumerateConfigurationsOptions const &options = *options_parser_ptr->value;

    auto make_enumerator_f = [&](std::string name, ConfigEnumInput const & initial_state) {
      return ConfigEnumRandomOccupations {
        initial_state,
        random_occ_params};
    };

    enumerate_configurations(
      options,
      make_enumerator_f,
      input_name_value_pairs.begin(),
      input_name_value_pairs.end(),
      primclex.db<Supercell>(),
      primclex.db<Configuration>(),
      logging);
  }

}
