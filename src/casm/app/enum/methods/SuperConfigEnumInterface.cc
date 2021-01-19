#include "casm/app/enum/methods/SuperConfigEnumInterface.hh"

#include "casm/app/APICommand.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler_impl.hh"
#include "casm/app/enum/dataformatter/ConfigEnumIO_impl.hh"
#include "casm/app/enum/enumerate_configurations_impl.hh"
#include "casm/app/enum/io/enumerate_configurations_json_io.hh"
#include "casm/app/enum/io/stream_io_impl.hh"
#include "casm/app/enum/standard_ConfigEnumInput_help.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/FillSupercell_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/SuperConfigEnum.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/io/SuperlatticeEnumeratorIO.hh"
#include "casm/database/io/json_io_impl.hh"

namespace CASM {

std::string SuperConfigEnumInterface::desc() const {
  std::string description =
      "The SuperConfigEnum method generates tilings of sub-configurations in "
      "selected \n"
      "supercells. The method works as follows: \n"
      "- Select a unit cell which all configurations tile into \n"
      "- Generate all unit-cell-filling equivalents of the input "
      "configurations. \n"
      "- Generate supercells of the unit cell and tile the supercells with all "
      "\n"
      "  combinations of the unit-cell-filling sub-configurations. Any global "
      "DoF are \n"
      "  set to zero value. \n\n";

  std::string custom_options =
      "  supercells: object, ScelEnum JSON settings (required)\n"
      "    Supercells to be tiled with selected configuration, specified in "
      "terms of\n"
      "    size and unit cell via a JSON object conforming to the format of "
      "'ScelEnum'\n"
      "    JSON settings \"min\", \"max\", \"dirs\", and \"unit_cell\". See "
      "'ScelEnum'\n"
      "    description for more details. The \"unit_cell\" must be a supercell "
      "of all   \n"
      "    all sub-configurations. \n\n"

      "  confignames: Array of strings (required) \n"
      "    Names of the sub-configurations that can be tiled into the "
      "\"unit_cell\".\n"
      "    Ex: \"confignames\" : "
      "[\"SCEL1_1_1_1_0_0_0/1\",\"SCEL2_2_1_1_0_0_0/3\"]\n\n"

      "  config_selection: string (optional) \n"
      "    Name of a selection of configurations to be tiled into selected "
      "supercells.\n\n";

  std::string examples =

      "  Examples:\n"
      "    To enumerate super-configurations of listed sub-configurations:\n"
      "     casm enum --method SuperConfigEnum -i \n"
      "     '{ \n"
      "        \"supercells\": { \n"
      "          \"max\": 4, \n"
      "          \"unit_cell\": \"SCEL2_1_2_1_0_0_0\" \n"
      "        },\n"
      "        \"confignames\": [\n"
      "          \"SCEL1_1_1_1_0_0_0/0\",\n"
      "          \"SCEL2_1_2_1_0_0_0/0\",\n"
      "          \"SCEL2_1_2_1_0_0_0/1\"\n"
      "        ]\n"
      "      }' \n"
      "\n"
      "    To enumerate super-configurations of listed sub-configurations from "
      "a \n"
      "    selection file:\n"
      "     casm enum --method SuperConfigEnum -i \n"
      "     '{ \n"
      "        \"supercells\": { \n"
      "          \"max\": 4, \n"
      "          \"unit_cell\": \"SCEL2_1_2_1_0_0_0\" \n"
      "        }, \n"
      "        \"config_selection\": \"selection_filename\"\n"
      "      }' \n"
      "\n";

  return name() + ": \n\n" + description + custom_options + examples;
}

std::string SuperConfigEnumInterface::name() const {
  return SuperConfigEnum::enumerator_name;
}

// check that the input configurations fit in "shared_unit_supercell"
void require_valid_sub_configurations(
    ParentInputParser &parser,
    DB::Selection<Configuration> const &config_selection,
    Supercell const &supercell) {
  auto it = config_selection.selected().begin();
  auto end = config_selection.selected().end();
  for (; it != end; ++it) {
    if (!is_valid_sub_configuration(it->ideal_lattice(), supercell)) {
      std::stringstream msg;
      msg << it.name() << " is not a valid sub configuration.";
      parser.error.insert(msg.str());
    }
  }
}

void SuperConfigEnumInterface::run(
    PrimClex &primclex, jsonParser const &json_options,
    jsonParser const &cli_options_as_json) const {
  Log &log = CASM::log();

  log.subsection().begin("SuperConfigEnum");
  ParentInputParser parser =
      make_enum_parent_parser(log, json_options, cli_options_as_json);
  std::runtime_error error_if_invalid{
      "Error reading SuperConfigEnum JSON input"};

  log.custom("Checking input");

  // 1) Parse ConfigEnumOptions ------------------
  auto options_parser_ptr = parser.parse_as<ConfigEnumOptions>(
      SuperConfigEnum::enumerator_name, primclex,
      primclex.settings().query_handler<Configuration>().dict());
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  ConfigEnumOptions const &options = *options_parser_ptr->value;
  print_options(log, options);
  log.set_verbosity(options.verbosity);

  // 2) Parse supercells ------------------

  auto shared_prim = primclex.shared_prim();
  auto scel_enum_props_subparser =
      parser.subparse_if<xtal::ScelEnumProps>("supercells");
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  // these parameters define which "target supercells" will be filled by
  // sub-configurations
  xtal::ScelEnumProps const &scel_enum_props =
      *scel_enum_props_subparser->value;

  // this is the smallest of the target supercells, the rest are supercells of
  // this
  auto shared_unit_supercell = std::make_shared<Supercell>(
      shared_prim, scel_enum_props.generating_matrix().cast<long>());

  // make target supercells (to be filled by sub-configurations to create
  // super-configurations)
  typedef std::string SupercellName;
  std::map<SupercellName, std::shared_ptr<Supercell const>> target_supercells;
  {
    ScelEnumByProps enumerator{shared_prim, scel_enum_props};
    for (auto const &supercell : enumerator) {
      target_supercells.emplace(supercell.name(),
                                std::make_shared<Supercell const>(supercell));
    }
  }
  print_initial_states(log, target_supercells);

  // 3) Parse sub-configurations --------------

  DB::Selection<Configuration> config_selection;
  try {
    config_selection = DB::make_selection<Configuration>(
        primclex.db<Configuration>(), parser.self, "confignames",
        "config_selection");
  } catch (std::exception &e) {
    std::stringstream msg;
    msg << "Error reading configurations: " << e.what();
    parser.error.insert(msg.str());
  }

  // logic check: require input configurations can fill the unit supercell
  require_valid_sub_configurations(parser, config_selection,
                                   *shared_unit_supercell);

  // if valid at this point, everything should be valid to make sub_configs and
  // target_supercells
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  // create the sub-configurations that will be building blocks of the
  // super-configurations note: the sub-configurations are themselves
  // super-configurations (of the input
  //       configurations) that fill the shared_unit_supercell
  std::vector<Configuration> sub_configs;
  {
    auto it = config_selection.selected().begin();
    auto end = config_selection.selected().end();
    Index n_before = 0;
    for (; it != end; ++it) {
      log << "Making equivalents of " << it.name() << ": " << std::flush;
      make_all_super_configurations(*it, shared_unit_supercell,
                                    std::back_inserter(sub_configs));
      log << sub_configs.size() - n_before << std::endl;
      n_before = sub_configs.size();
    }
    // TODO: check for and remove duplicates caused by non-primitive input
    // configurations?
    log << "Total sub-configurations, including all equivalents: "
        << sub_configs.size() << std::flush;
  }

  // 4) Enumerate configurations ------------------

  // Functor to construct SuperConfigEnum for each target supercell
  auto make_enumerator_f =
      [&](Index index, std::string name,
          std::shared_ptr<Supercell const> const &target_supercell) {
        return SuperConfigEnum{*target_supercell, sub_configs.begin(),
                               sub_configs.end()};
      };

  typedef ConfigEnumData<SuperConfigEnum, std::shared_ptr<Supercell const>>
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
  formatter.push_back(ConfigEnumIO::initial_state_index<ConfigEnumDataType>(),
                      ConfigEnumIO::initial_state_name<ConfigEnumDataType>());

  log << std::endl;
  log.begin("SuperConfigEnum enumeration");

  enumerate_configurations(primclex, options, make_enumerator_f,
                           target_supercells.begin(), target_supercells.end(),
                           formatter);

  log.end_section();
}

}  // namespace CASM
