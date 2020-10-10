#include "casm/app/enum/enumerate_configurations_impl.hh"
#include "casm/app/enum/methods/ConfigEnumAllOccupationsInterface.hh"
#include "casm/clex/ConfigEnumAllOccupations.hh"

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

  std::string ConfigEnumAllOccupationsInterface::desc() const {

    std::string custom_options = "";

    std::string examples =
      "  Examples:\n"
      "    To enumerate all occupations in supercells up to and including size 4:\n"
      "      casm enum --method ConfigEnumAllOccupations -i '{\"supercells\": {\"max\": 4}}' \n"
      "\n"
      "    To enumerate all occupations in all existing supercells:\n"
      "      casm enum --method ConfigEnumAllOccupations --all\n"
      "\n"
      "    To enumerate all occupations in particular supercells:\n"
      "      casm enum --method ConfigEnumAllOccupations -i \n"
      "      '{ \n"
      "        \"scelnames\": [\n"
      "           \"SCEL1_1_1_1_0_0_0\",\n"
      "           \"SCEL2_1_2_1_0_0_0\",\n"
      "           \"SCEL4_1_4_1_0_0_0\"\n"
      "        ]\n"
      "      }' \n\n";

    return name() + ": \n\n" + custom_options + standard_ConfigEnumInput_help() + examples;
  }

  std::string ConfigEnumAllOccupationsInterface::name() const {
    return ConfigEnumAllOccupations::enumerator_name;
  }

  void ConfigEnumAllOccupationsInterface::run(
    PrimClex &primclex,
    jsonParser const &json_options,
    jsonParser const &cli_options_as_json) const {

    Log &log = CASM::log();

    // combine JSON options and CLI options
    jsonParser json_combined = combine_configuration_enum_json_options(
                                 json_options,
                                 cli_options_as_json);

    // Read input data from JSON
    ParentInputParser parser {json_combined};

    auto input_parser_ptr = parser.parse_as<std::vector<std::pair<std::string, ConfigEnumInput>>>(
                              primclex.shared_prim(),
                              &primclex,
                              primclex.db<Supercell>(),
                              primclex.db<Configuration>());
    auto options_parser_ptr = parser.parse_as<EnumerateConfigurationsOptions>(
                                ConfigEnumAllOccupations::enumerator_name,
                                primclex,
                                primclex.settings().query_handler<Configuration>().dict());

    std::runtime_error error_if_invalid {"Error reading ConfigEnumAllOccupations JSON input"};
    report_and_throw_if_invalid(parser, log, error_if_invalid);

    auto const &input_name_value_pairs = *input_parser_ptr->value;
    EnumerateConfigurationsOptions const &options = *options_parser_ptr->value;

    auto make_enumerator_f = [&](std::string name, ConfigEnumInput const & initial_state) {
      return ConfigEnumAllOccupations {initial_state};
    };

    enumerate_configurations(
      options,
      make_enumerator_f,
      input_name_value_pairs.begin(),
      input_name_value_pairs.end(),
      primclex.db<Supercell>(),
      primclex.db<Configuration>());
  }

}
