#include "casm/app/enum/methods/ConfigEnumInterfaceTemplate.hh"

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
#include "casm/clex/PrimClex.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/InputEnumerator.hh"
#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"

/// To add a new `casm enum` method:
/// Change "ConfigEnumTypeTemplate" to "ConfigEnum<Something>"
/// Change "ConfigEnumInterfaceTemplate" to "ConfigEnum<Something>Interface"
/// Look for the **EDIT THIS** tags for customizing the implementation

namespace CASM {

/// This class implements a method to enumerate Configuration
/// If useful outside of `casm enum`, it should be placed in casm/clex.
class ConfigEnumTypeTemplate : public InputEnumeratorBase<Configuration> {
  // -- Required members -------------------

 public:
  ConfigEnumTypeTemplate(ConfigEnumInput const &_in_config, Index _n_config,
                         MTRand &_mtrand);

  std::string name() const override;

  static const std::string enumerator_name;

 private:
  /// Implements increment
  void increment() override;

  // **EDIT THIS**, the following custom methods and members
  // should be edited to implement the enumeration method
  // -- Custom members -------------------

  void randomize() {
    // for each selected site, choose a random occupant
    for (Index l : m_site_selection) {
      m_current->set_occ(l, m_mtrand.randInt(m_max_allowed[l]));
    }
  }

  Index m_n_config;
  MTRand &m_mtrand;
  Eigen::VectorXi m_max_allowed;
  std::vector<Index> m_site_selection;
  notstd::cloneable_ptr<Configuration> m_current;
};

// **EDIT THIS**, to set up the enumeration method
/// \brief Constructor
///
/// The constructor must:
/// 1) Set up the enumeration method
/// 2) Call this->_initialize(Configuration const *) to initialize
/// the enumeration and point at the Configuration that will be
/// returned when dereferencing the `begin` iterator.
/// or 2) call this->_invalidate() if the input parameters are such
/// that no Configuration can be enumerated
///
ConfigEnumTypeTemplate::ConfigEnumTypeTemplate(
    ConfigEnumInput const &_in_config, Index _n_config, MTRand &_mtrand)
    : m_n_config(_n_config),
      m_mtrand(_mtrand),
      m_max_allowed(
          _in_config.configuration().supercell().max_allowed_occupation()),
      m_site_selection(_in_config.sites().begin(), _in_config.sites().end()) {
  if (m_n_config < 0) {
    throw std::runtime_error("Error in ConfigEnumTypeTemplate: n_config < 0");
  }
  if (m_n_config == 0) {
    this->_invalidate();
    return;
  }

  // Copy the input configuration
  m_current = notstd::make_cloneable<Configuration>(_in_config.configuration());

  // Make sure it has no properties
  reset_properties(*m_current);

  this->_initialize(&(*m_current));

  // Make initial random config
  this->randomize();
}

std::string ConfigEnumTypeTemplate::name() const { return enumerator_name; }

// **EDIT THIS**, to return the correct name
const std::string ConfigEnumTypeTemplate::enumerator_name =
    "ConfigEnumTypeTemplate";

// **EDIT THIS**, to generate the next Configuration
/// The increment method must:
/// 1) Call this->_increment_step()
/// 2) Do something to change m_current
/// or 2) Call this->_invalidate() if there are no more Configurations
/// to enumerate by this method
///
/// Note: this->step() gives the current enumeration step
/// Note: if dereferencing an iterator should return a reference to a
/// different configuration than m_current, call
/// this->_set_current_ptr(Configuration const *) to specify which
/// Configuration should be referenced
void ConfigEnumTypeTemplate::increment() {
  this->_increment_step();
  if (step() < m_n_config) {
    // Make next random config
    this->randomize();
  } else {
    this->_invalidate();
  }
}

std::string ConfigEnumInterfaceTemplate::desc() const {
  // **EDIT THIS**, to update the description of the method
  std::string description = "Generate random occupation orderings.\n\n";

  // **EDIT THIS**, if there are custom options
  std::string custom_options =
      "  n_config: integer (required) \n"
      "    How many random configurations to generate. Includes duplicate\n"
      "    and pre- existing configurations. \n\n"

      "  option_a: integer (optional, default=0) \n"
      "    An optional parameter, with default value 0. \n\n"

      "  option_b: integer (optional) \n"
      "    Another optional parameter. \n\n";

  // **EDIT THIS**, to provide the user some example input
  std::string examples =
      "  Examples:\n"
      "    To enumerate 200 random occupations in supercells up to and "
      "including size 4:\n"
      "      casm enum --method ConfigEnumRandomOccupations -i \n"
      "        '{\"supercells\":{\"max\":4}, \"n_config\": 200}' \n"
      "\n"
      "    To enumerate 200 random occupations in all existing supercells:\n"
      "      casm enum --method ConfigEnumRandomOccupations --all -i "
      "'{\"n_config\": 200}' \n"
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

  return name() + ": \n\n" + description + custom_options +
         standard_ConfigEnumInput_help() + examples;
}

std::string ConfigEnumInterfaceTemplate::name() const {
  return ConfigEnumTypeTemplate::enumerator_name;
}

/// This parses input parameters and runs the enumeration method
void ConfigEnumInterfaceTemplate::run(
    PrimClex &primclex, jsonParser const &json_options,
    jsonParser const &cli_options_as_json) const {
  Log &log = CASM::log();

  /// This combines JSON input (from -i or -s) with JSON
  /// representing other cli input options
  log.subsection().begin("ConfigEnumInterfaceTemplate");
  ParentInputParser parser =
      make_enum_parent_parser(log, json_options, cli_options_as_json);
  std::runtime_error error_if_invalid{
      "Error reading ConfigEnumInterfaceTemplate JSON input"};

  log.custom("Checking input");

  // 1) Parse ConfigEnumOptions ------------------

  /// This parses the standard options that get stored as
  /// ConfigEnumOptions. These options get used automatically
  /// in `enumerate_configurations` and for generating
  /// user-specified data formatters
  auto options_parser_ptr = parser.parse_as<ConfigEnumOptions>(
      ConfigEnumTypeTemplate::enumerator_name, primclex,
      primclex.settings().query_handler<Configuration>().dict());
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  ConfigEnumOptions const &options = *options_parser_ptr->value;
  print_options(log, options);
  log.set_verbosity(options.verbosity);

  // 2) Parse initial enumeration states ------------------

  /// This reads the standard input options for specifying the initial
  /// enumeration states (ConfigEnumInput, a combination of
  /// configuration & selected sites).
  auto input_parser_ptr =
      parser.parse_as<std::vector<std::pair<std::string, ConfigEnumInput>>>(
          primclex.shared_prim(), &primclex, primclex.db<Supercell>(),
          primclex.db<Configuration>());
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  auto const &named_initial_states = *input_parser_ptr->value;
  print_initial_states(log, named_initial_states);

  //**EDIT THIS**, edit this section to read method parameters from JSON
  // 3) Parse custom parameters ------------------

  // random number generator, default initialization
  MTRand mtrand;

  // example parsing a required input,
  // the number of random configuration to generate
  Index n_config;
  parser.require(n_config, "n_config");

  // example error message for input parameter
  if (n_config <= 0) {
    fs::path option = "n_config";
    std::string message = "Error: n_config must be > 0";
    parser.insert_error(option, message);
  }

  // example parsing an optional input, with default value 0
  // no error message if "option_a" is not included in the input,
  // but an error message is provided if input value
  // cannot be converted to Index
  Index an_optional_parameter = 0;
  parser.optional(an_optional_parameter, "option_a");

  // another example parsing an optional input, with default nullopt
  // no error message if "option_b" is not included in the input,
  // but an error message is provided if input value
  // cannot be converted to Index
  std::optional<Index> another_optional_parameter = std::nullopt;
  parser.optional(another_optional_parameter, "option_b");

  // example warning message for input parameter
  if (n_config > 1e6) {
    fs::path option = "n_config";
    std::string message =
        "Error: n_config is very large, this may take a long time";
    parser.insert_warning(option, message);
  }

  // this will throw a nice error message if the required options are
  // not given, if values are the wrong type, etc.
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  // 4) Enumerate configurations ------------------

  // **EDIT THIS** to pass the input options and construct the
  // ConfigEnumTypeTemplate enumeration method.

  // This is a functor that creates a new enumerator for each initial state
  // (ConfigEnumInput, a combination of configuration & selected sites)
  // specified by the input file.
  auto make_enumerator_f = [&](Index index, std::string name,
                               ConfigEnumInput const &initial_state) {
    return ConfigEnumTypeTemplate{initial_state, n_config, mtrand};
  };

  // Create a DataFormatter<Configuration> to output data about
  // the enumerated configurations if "output_configurations"==true
  typedef ConfigEnumData<ConfigEnumTypeTemplate, ConfigEnumInput>
      ConfigEnumDataType;
  DataFormatter<ConfigEnumDataType> formatter;

  // These formatters should be appropriate for most Configuration enumerators
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

  // This adds `casm query` property formatters specified by the user
  // via "output_configurations_options"/"properties"
  for (const auto &formatter_ptr : options.output_formatter.formatters()) {
    formatter.push_back(
        make_datum_formatter_adapter<ConfigEnumDataType, Configuration>(
            *formatter_ptr));
  }

  log << std::endl;
  log.begin("ConfigEnumTypeTemplate enumeration");

  // This runs the enumeration for all initial states, stores
  // and output results
  enumerate_configurations(primclex, options, make_enumerator_f,
                           named_initial_states.begin(),
                           named_initial_states.end(), formatter);

  log.end_section();
}

}  // namespace CASM
