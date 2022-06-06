#include "casm/app/enum/methods/ConfigEnumRandomLocalInterface.hh"

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
#include "casm/clex/ConfigEnumRandomLocal.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"

namespace CASM {

std::string ConfigEnumRandomLocalInterface::desc() const {
  std::string custom_options =
      "  n_config: integer (optional, default=100) \n"
      "    How many random configurations to generate. Includes duplicate and "
      "pre-\n"
      "    existing configurations.                                            "
      "     \n\n"

      "  magnitude: positive number (optional, default=1.0) \n"
      "    Magnitude used to scale random vector at each site. If "
      "\"distribution\" == \"normal\",\n"
      "    magnitude specifies standard deviation of D-dimensional Gaussian (D "
      "is dimension of \n"
      "    site DoF value). If \"distribution\" == \"uniform\", magnitude is "
      "radius of D-dimensional\n"
      "    ball from which random vectors are chosen.\n\n"

      "  dof: string (required) \n"
      "    Name of site degree of freecom  for which normal coordinates are to "
      "be generated.\n"
      "    Must be one of the degrees of freedom under consideration in the "
      "current project,\n"
      "    as determined by prim.json\n\n"

      "  distribution: string (optional, default=\"normal\") \n"
      "    Distribution from which perturbation vectors are drawn. Options are "
      "\"uniform\"\n"
      "    (i.e., uniform on the unit sphere), or \"normal\" (i.e., "
      "zero-centered Gaussian).\n\n";

  std::string examples =
      "  Examples:\n"
      "    To enumerate 200 random occupations in supercells up to and "
      "including size 4:\n"
      "      casm enum --method ConfigEnumRandomLocal -i \n"
      "        '{\"supercells\":{\"max\":4}, \"n_config\": 200}' \n"
      "\n"
      "    To enumerate 200 random occupations in all existing supercells:\n"
      "      casm enum --method ConfigEnumRandomLocal -i '{\"n_config\": 200}' "
      "\n"
      "\n"
      "    To enumerate 100 random occupations in all existing supercells:\n"
      "      casm enum --method ConfigEnumRandomLocal --all \n"
      "\n"
      "    To enumerate 200 random occupations in particular supercells:\n"
      "      casm enum --method ConfigEnumRandomLocal -i \n"
      "      '{ \n"
      "        \"scelnames\": [\n"
      "           \"SCEL1_1_1_1_0_0_0\",\n"
      "           \"SCEL2_1_2_1_0_0_0\",\n"
      "           \"SCEL4_1_4_1_0_0_0\"\n"
      "        ],\n"
      "        \"n_config\": 200\n"
      "      }' \n\n";

  return name() + ": \n\n" + custom_options + standard_ConfigEnumInput_help() +
         examples;
}

std::string ConfigEnumRandomLocalInterface::name() const {
  return ConfigEnumRandomLocal::enumerator_name;
}

void parse(InputParser<ConfigEnumRandomLocalParams> &parser, MTRand &mtrand) {
  DoFKey dof_key;
  double mag;
  Index n_config;
  bool normal_distribution;

  parser.require(dof_key, "dof");
  parser.optional_else(mag, "magnitude", double{1.0});
  parser.optional_else(n_config, "n_config", Index{100});

  std::string distribution_name;
  parser.optional_else(distribution_name, "distribution",
                       std::string{"normal"});
  if (distribution_name == "normal") {
    normal_distribution = true;
  } else if (distribution_name == "uniform") {
    normal_distribution = false;
  } else {
    std::stringstream msg;
    msg << "Error: \"distribution\" must be \"normal\" (default) or "
           "\"uniform\"";
    parser.error.insert(msg.str());
  }

  parser.value = notstd::make_unique<ConfigEnumRandomLocalParams>(
      mtrand, dof_key, n_config, mag, normal_distribution);
}

void ConfigEnumRandomLocalInterface::run(
    PrimClex &primclex, jsonParser const &json_options,
    jsonParser const &cli_options_as_json) const {
  Log &log = CASM::log();

  log.subsection().begin("ConfigEnumRandomLocal");
  ParentInputParser parser =
      make_enum_parent_parser(log, json_options, cli_options_as_json);
  std::runtime_error error_if_invalid{
      "Error reading ConfigEnumRandomLocal JSON input"};

  log.custom("Checking input");

  // 1) Parse ConfigEnumOptions ------------------

  auto options_parser_ptr = parser.parse_as<ConfigEnumOptions>(
      ConfigEnumRandomLocal::enumerator_name, primclex,
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

  // 3) Parse ConfigEnumRandomLocalParams ------------------

  MTRand mtrand;
  auto random_local_params_parser_ptr =
      parser.parse_as<ConfigEnumRandomLocalParams>(mtrand);
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  ConfigEnumRandomLocalParams const &random_local_params =
      *random_local_params_parser_ptr->value;

  // 4) Enumerate configurations ------------------

  auto make_enumerator_f = [&](Index index, std::string name,
                               ConfigEnumInput const &initial_state) {
    return ConfigEnumRandomLocal{initial_state, random_local_params};
  };

  typedef ConfigEnumData<ConfigEnumRandomLocal, ConfigEnumInput>
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
  log.begin("ConfigEnumRandomLocal enumeration");

  enumerate_configurations(primclex, options, make_enumerator_f,
                           named_initial_states.begin(),
                           named_initial_states.end(), formatter);

  log.end_section();
}

}  // namespace CASM
