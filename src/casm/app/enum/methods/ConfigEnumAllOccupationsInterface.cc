#include "casm/app/enum/methods/ConfigEnumAllOccupationsInterface.hh"

#include "casm/app/APICommand.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler_impl.hh"
#include "casm/app/enum/dataformatter/ConfigEnumIO_impl.hh"
#include "casm/app/enum/enumerate_configurations_impl.hh"
#include "casm/app/enum/io/enumerate_configurations_json_io.hh"
#include "casm/app/enum/io/stream_io_impl.hh"
#include "casm/app/enum/standard_ConfigEnumInput_help.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"

namespace CASM {

std::string ConfigEnumAllOccupationsInterface::desc() const {
  std::string custom_options = "";

  std::string examples =
      "  Examples:\n"
      "    To enumerate all occupations in supercells up to and including size "
      "4:\n"
      "      casm enum --method ConfigEnumAllOccupations -i '{\"supercells\": "
      "{\"max\": 4}}' \n"
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

  return name() + ": \n\n" + custom_options + standard_ConfigEnumInput_help() +
         examples;
}

std::string ConfigEnumAllOccupationsInterface::name() const {
  return ConfigEnumAllOccupations::enumerator_name;
}

void ConfigEnumAllOccupationsInterface::run(
    PrimClex &primclex, jsonParser const &json_options,
    jsonParser const &cli_options_as_json) const {
  Log &log = CASM::log();

  log.subsection().begin("ConfigEnumAllOccupations");
  ParentInputParser parser =
      make_enum_parent_parser(log, json_options, cli_options_as_json);
  std::runtime_error error_if_invalid{
      "Error reading ConfigEnumAllOccupations JSON input"};

  log.custom("Checking input");

  // 1) Parse ConfigEnumOptions ------------------

  auto options_parser_ptr = parser.parse_as<ConfigEnumOptions>(
      ConfigEnumAllOccupations::enumerator_name, primclex,
      primclex.settings().query_handler<Configuration>().dict());
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  ConfigEnumOptions const &options = *options_parser_ptr->value;
  print_options(log, options);
  log.set_verbosity(options.verbosity);

  // 2) Parse initial enumeration states ------------------

  auto input_parser_ptr =
      parser.parse_as<std::vector<std::pair<std::string, ConfigEnumInput>>>(
          primclex.shared_prim(), &primclex, primclex.db<Supercell>(),
          primclex.db<Configuration>());
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  auto const &named_initial_states = *input_parser_ptr->value;
  print_initial_states(log, named_initial_states);

  // 3) Enumerate configurations ------------------

  auto make_enumerator_f = [&](Index index, std::string name,
                               ConfigEnumInput const &initial_state) {
    return ConfigEnumAllOccupations{initial_state};
  };

  typedef ConfigEnumData<ConfigEnumAllOccupations, ConfigEnumInput>
      ConfigEnumDataType;
  DataFormatter<ConfigEnumDataType> formatter;
  formatter.push_back(ConfigEnumIO::name<ConfigEnumDataType>(),
                      ConfigEnumIO::selected<ConfigEnumDataType>(),
                      ConfigEnumIO::is_new<ConfigEnumDataType>(),
                      ConfigEnumIO::is_existing<ConfigEnumDataType>());
  if (options.filter) {
    formatter.push_back(
        ConfigEnumIO::is_excluded_by_filter<ConfigEnumDataType>());
  }
  formatter.push_back(
      ConfigEnumIO::initial_state_index<ConfigEnumDataType>(),
      ConfigEnumIO::initial_state_name<ConfigEnumDataType>(),
      ConfigEnumIO::initial_state_configname<ConfigEnumDataType>(),
      ConfigEnumIO::n_selected_sites<ConfigEnumDataType>());

  log << std::endl;
  log.begin("ConfigEnumAllOccupations enumeration");

  enumerate_configurations(primclex, options, make_enumerator_f,
                           named_initial_states.begin(),
                           named_initial_states.end(), formatter);

  log.end_section();
}

}  // namespace CASM
