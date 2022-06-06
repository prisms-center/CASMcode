#include "casm/app/enum/methods/ConfigEnumAllOccupationsInterface.hh"

#include "casm/app/APICommand.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler_impl.hh"
#include "casm/app/enum/dataformatter/ConfigEnumIO_impl.hh"
#include "casm/app/enum/enumerate_configurations_impl.hh"
#include "casm/app/enum/io/enumerate_configurations_json_io.hh"
#include "casm/app/enum/io/stream_io_impl.hh"
#include "casm/app/enum/standard_ConfigEnumInput_help.hh"
#include "casm/casm_io/dataformatter/DatumFormatterAdapter.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/optional.hh"
#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"

namespace CASM {

std::string ConfigEnumAllOccupationsInterface::desc() const {
  std::string custom_options =
      "  skip_non_primitive: bool (optional)\n"
      "    If true, non-primitive configurations are skipped in  \n"
      "    the enumeration. If false, they are included. Whether \n"
      "    they are included in the database or not is specified \n"
      "    separately by the \"primitive_only\" option. This     \n"
      "    option allows including non-primitive configurations  \n"
      "    in the output generated when \"output_configurations\"==true.\n"
      "    The default value is true if enumeration is occuring  \n"
      "    on all sites in the configuration, and false otherwise.\n\n"

      "  skip_non_canonical: bool (optional)\n"
      "    If true, non-canonical configurations are skipped in  \n"
      "    the enumeration. If false, they are included. For     \n"
      "    either choice, only canonical configurations are ever \n"
      "    inserted into the configuration database. This option \n"
      "    allows including non-canonical configurations in the  \n"
      "    output generated when \"output_configurations\"==true.\n"
      "    The default value is true if enumeration is occuring  \n"
      "    on all sites in the configuration, and false otherwise.\n\n";

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

  return name() + ": \n\n" + standard_ConfigEnumInput_help() + custom_options +
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

  // 1a) Parse ConfigEnumOptions ------------------

  auto options_parser_ptr = parser.parse_as<ConfigEnumOptions>(
      ConfigEnumAllOccupations::enumerator_name, primclex,
      primclex.settings().query_handler<Configuration>().dict());
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  ConfigEnumOptions const &options = *options_parser_ptr->value;
  print_options(log, options);
  log.set_verbosity(options.verbosity);

  // 1b) Parse custom options ---------------------
  std::optional<bool> skip_non_primitive;
  parser.optional(skip_non_primitive, "skip_non_primitive");

  std::optional<bool> skip_non_canonical;
  parser.optional(skip_non_canonical, "skip_non_canonical");

  log << std::boolalpha;
  log.indent() << "skip_non_primitive:";
  if (skip_non_primitive.has_value()) {
    log << skip_non_primitive.value() << std::endl;
  } else {
    log << "null" << std::endl;
  }
  log.indent() << "skip_non_canonical:";
  if (skip_non_canonical.has_value()) {
    log << skip_non_canonical.value() << std::endl;
  } else {
    log << "null" << std::endl;
  }
  log << std::noboolalpha;

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
    // set defaults
    bool primitive_only = true;
    bool canonical_only = true;
    if (initial_state.sites().size() != initial_state.configuration().size()) {
      primitive_only = false;
      canonical_only = false;
    }

    // optionally, override defaults
    if (skip_non_primitive.has_value()) {
      primitive_only = skip_non_primitive.value();
    }
    if (skip_non_canonical.has_value()) {
      canonical_only = skip_non_canonical.value();
    }
    return ConfigEnumAllOccupations{initial_state, primitive_only,
                                    canonical_only};
  };

  typedef ConfigEnumData<ConfigEnumAllOccupations, ConfigEnumInput>
      ConfigEnumDataType;
  DataFormatter<ConfigEnumDataType> formatter;
  formatter.push_back(ConfigEnumIO::canonical_configname<ConfigEnumDataType>(),
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
  for (const auto &formatter_ptr : options.output_formatter.formatters()) {
    formatter.push_back(
        make_datum_formatter_adapter<ConfigEnumDataType, Configuration>(
            *formatter_ptr));
  }

  log << std::endl;
  log.begin("ConfigEnumAllOccupations enumeration");

  enumerate_configurations(primclex, options, make_enumerator_f,
                           named_initial_states.begin(),
                           named_initial_states.end(), formatter);

  log.end_section();
}

}  // namespace CASM
