// This is an example enumerator plugin implementation. It is a copy of ConfigEnumAllOccupations,
// with the name changed to TestEnum. For a plugin, the primary requirements are:
// - Must have a source file named "MethodName.cc" in the enumerator plugins directory
// - Must have an interface generating function. It should be:
//   - An `extern "C"` function
//   - With name and signature `CASM::EnumInterfaceBase *make_<MethodName>_interface()`
// - Placing declarations in a header file is optional
// - EnumInterfaceBase is defined in "casm/app/enum/EnumInterface.hh"
// - See the standard interfaces implemented in `casm/app/enum/methods` for examples

// --- TestEnum ---

#include "casm/casm_io/Log.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/container/Counter.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/InputEnumerator.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

  /// Enumerate over all possible occupations on particular sites in a Configuration
  ///
  class TestEnum : public InputEnumeratorBase<Configuration> {

    // -- Required members -------------------

  public:

    /// Construct with a ConfigEnumInput, specifying which sites to enumerate on and which to keep fixed
    ///
    /// Note:
    /// - The output configurations are set to fixed occupation matching
    ///   `config_enum_input.configuration()` on all sites not included in
    ///   `config_enum_input.sites()`.
    /// - All allowed occupation values are enumerated on the selected sites
    ///   (`config_enum_input.sites()`), but Configuration are only output if the are primitive.
    /// - If all sites are selected, then only canonical Configuration are output. This can be
    ///   checked with `this->canonical_guarantee()`.
    TestEnum(ConfigEnumInput const &config_enum_input);

    std::string name() const override;

    /// Returns true if enumerator is guaranteed to output canonical configurations
    bool canonical_guarantee() const;

    static const std::string enumerator_name;

  private:

    /// Implements increment
    void increment() override;

    // -- Unique -------------------

    /// Returns true if current() is valid for output
    ///
    /// When enumerating on all sites:
    /// - Require output configurations are primitive and canonical
    /// When enumerating on a subset of sites:
    /// - Require output configurations are primitive (but not canonical)
    ///
    bool _current_is_valid_for_output() const;

    /// Site index to enumerate on
    std::set<Index> m_site_index_selection;

    /// Counter over allowed occupation indices on sites in m_site_index_selection
    Counter<std::vector<int> > m_counter;

    /// The current configuration
    notstd::cloneable_ptr<Configuration> m_current;

    /// True if enumerating on a subset of supercell sites
    bool m_enumerate_on_a_subset_of_supercell_sites;
  };

  namespace local_impl {
    std::vector<int> max_selected_occupation(ConfigEnumInput const &config_enum_input) {
      auto const &supercell = config_enum_input.configuration().supercell();
      std::vector<int> max_allowed = supercell.max_allowed_occupation();

      std::vector<int> max_allowed_on_selected_sites;
      for(Index i : config_enum_input.sites()) {
        max_allowed_on_selected_sites.push_back(max_allowed[i]);
      }

      return max_allowed_on_selected_sites;
    }
  }

  /// Conditionally true for TestEnum (true when enumerating on all sites)
  template<>
  bool is_guaranteed_for_database_insert(TestEnum const &enumerator) {
    return enumerator.canonical_guarantee();
  }

  /// \brief Construct with a Supercell, using all permutations
  TestEnum::TestEnum(const ConfigEnumInput &config_enum_input) :
    m_site_index_selection(config_enum_input.sites()),
    m_counter(
      std::vector<int>(config_enum_input.sites().size(), 0),
      local_impl::max_selected_occupation(config_enum_input),
      std::vector<int>(config_enum_input.sites().size(), 1)),
    m_current(notstd::make_cloneable<Configuration>(config_enum_input.configuration())),
    m_enumerate_on_a_subset_of_supercell_sites(
      m_site_index_selection.size() != config_enum_input.configuration().size()) {

    m_current->set_occupation(m_counter());
    reset_properties(*m_current);
    this->_initialize(&(*m_current));

    // Make sure that current() is a primitive canonical config
    if(!_current_is_valid_for_output()) {
      increment();
    }

    // set step to 0
    if(valid()) {
      _set_step(0);
    }
    m_current->set_source(this->source(step()));
  }

  std::string TestEnum::name() const {
    return enumerator_name;
  }

  /// Returns true if enumerator is guaranteed to output canonical configurations
  ///
  /// Output is sometimes canonical (when enumerating on all supercell sites) and sometimes not
  /// (when enumerating on a subset of supercell sites). If this returns true, then the output is
  /// guaranteed to be in canonical form.
  bool TestEnum::canonical_guarantee() const {
    return !m_enumerate_on_a_subset_of_supercell_sites;
  }

  const std::string TestEnum::enumerator_name = "TestEnum";

  /// Implements _increment over all occupations
  void TestEnum::increment() {

    bool is_valid_config {false};

    while(!is_valid_config && ++m_counter) {
      for(Index site_index : m_site_index_selection) {
        m_current->set_occ(site_index, m_counter[site_index]);
      }
      is_valid_config = _current_is_valid_for_output();
    }

    // if while loop exited with valid configuration, increment step number
    if(m_counter.valid()) {
      this->_increment_step();
    }
    // if while loop exited because counter hit end, invalidate this enumerator
    else {
      this->_invalidate();
    }
    m_current->set_source(this->source(step()));
  }

  /// Returns true if current() is primitive and canonical
  bool TestEnum::_current_is_valid_for_output() const {
    if(m_enumerate_on_a_subset_of_supercell_sites) {
      return current().is_primitive();
    }
    else {
      return current().is_primitive() && current().is_canonical();
    }
  }

}

// --- TestEnumInterface ---
#include "casm/app/enum/enumerate_configurations_impl.hh"

#include "casm/app/APICommand.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler_impl.hh"
#include "casm/app/enum/standard_ConfigEnumInput_help.hh"
#include "casm/app/enum/dataformatter/ConfigEnumIO_impl.hh"
#include "casm/app/enum/io/enumerate_configurations_json_io.hh"
#include "casm/app/enum/io/stream_io_impl.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"

namespace CASM {

  class TestEnumInterface : public EnumInterfaceBase {
    CLONEABLE(TestEnumInterface)
  public:

    std::string desc() const override {
      std::string custom_options = "";

      std::string examples =
        "  Examples:\n"
        "    To enumerate all occupations in supercells up to and including size 4:\n"
        "      casm enum --method TestEnum -i '{\"supercells\": {\"max\": 4}}' \n"
        "\n"
        "    To enumerate all occupations in all existing supercells:\n"
        "      casm enum --method TestEnum --all\n"
        "\n"
        "    To enumerate all occupations in particular supercells:\n"
        "      casm enum --method TestEnum -i \n"
        "      '{ \n"
        "        \"scelnames\": [\n"
        "           \"SCEL1_1_1_1_0_0_0\",\n"
        "           \"SCEL2_1_2_1_0_0_0\",\n"
        "           \"SCEL4_1_4_1_0_0_0\"\n"
        "        ]\n"
        "      }' \n\n";

      return name() + ": \n\n" + custom_options + standard_ConfigEnumInput_help() + examples;
    }

    std::string name() const override {
      return TestEnum::enumerator_name;
    }

    void run(PrimClex &primclex,
             jsonParser const &json_options,
             jsonParser const &cli_options_as_json) const override {

      Log &log = CASM::log();

      log.subsection().begin("TestEnum");
      ParentInputParser parser = make_enum_parent_parser(log, json_options, cli_options_as_json);
      std::runtime_error error_if_invalid {"Error reading TestEnum JSON input"};

      log.custom("Checking input");

      // 1) Parse ConfigEnumOptions ------------------

      auto options_parser_ptr = parser.parse_as<ConfigEnumOptions>(
                                  TestEnum::enumerator_name,
                                  primclex,
                                  primclex.settings().query_handler<Configuration>().dict());
      report_and_throw_if_invalid(parser, log, error_if_invalid);
      ConfigEnumOptions const &options = *options_parser_ptr->value;
      print_options(log, options);
      log.set_verbosity(options.verbosity);

      // 2) Parse initial enumeration states ------------------

      auto input_parser_ptr = parser.parse_as<std::vector<std::pair<std::string, ConfigEnumInput>>>(
                                primclex.shared_prim(),
                                &primclex,
                                primclex.db<Supercell>(),
                                primclex.db<Configuration>());
      report_and_throw_if_invalid(parser, log, error_if_invalid);
      auto const &named_initial_states = *input_parser_ptr->value;
      print_initial_states(log, named_initial_states);

      // 3) Enumerate configurations ------------------

      auto make_enumerator_f = [&](Index index, std::string name, ConfigEnumInput const & initial_state) {
        return TestEnum {initial_state};
      };

      typedef ConfigEnumData<TestEnum, ConfigEnumInput> ConfigEnumDataType;
      DataFormatter<ConfigEnumDataType> formatter;
      formatter.push_back(
        ConfigEnumIO::name<ConfigEnumDataType>(),
        ConfigEnumIO::selected<ConfigEnumDataType>(),
        ConfigEnumIO::is_new<ConfigEnumDataType>(),
        ConfigEnumIO::is_existing<ConfigEnumDataType>()
      );
      if(options.filter) {
        formatter.push_back(ConfigEnumIO::is_excluded_by_filter<ConfigEnumDataType>());
      }
      formatter.push_back(
        ConfigEnumIO::initial_state_index<ConfigEnumDataType>(),
        ConfigEnumIO::initial_state_name<ConfigEnumDataType>(),
        ConfigEnumIO::initial_state_configname<ConfigEnumDataType>(),
        ConfigEnumIO::n_selected_sites<ConfigEnumDataType>()
      );

      log << std::endl;
      log.begin("TestEnum enumeration");

      enumerate_configurations(
        primclex,
        options,
        make_enumerator_f,
        named_initial_states.begin(),
        named_initial_states.end(),
        formatter);

      log.end_section();
    }

  };
}

// --- make_TestEnum_interface ---

extern "C" {
  CASM::EnumInterfaceBase *make_TestEnum_interface() {
    return new CASM::TestEnumInterface {};
  }
}
