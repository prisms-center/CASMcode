#include "casm/app/enum/enumerate_configurations_impl.hh"
#include "casm/app/enum/methods/ConfigEnumStrainInterface.hh"
#include "casm/clex/ConfigEnumStrain.hh"

#include "casm/app/APICommand.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler_impl.hh"
#include "casm/app/enum/standard_ConfigEnumInput_help.hh"
#include "casm/app/enum/dataformatter/ConfigEnumIO_impl.hh"
#include "casm/app/enum/io/enumerate_configurations_json_io.hh"
#include "casm/app/enum/io/stream_io_impl.hh"
#include "casm/casm_io/dataformatter/DatumFormatterAdapter.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/ConfigIOStrain.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/DoFSpace.hh"
#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"
#include "casm/symmetry/SymRepTools.hh"
#include "casm/symmetry/io/json/SymRepTools.hh"

namespace CASM {

  // namespace adapter {
  //
  //   Configuration const &
  //   Adapter<Configuration, ConfigEnumData<ConfigEnumStrain, ConfigEnumInput>>::operator()(
  //     ConfigEnumData<ConfigEnumStrain, ConfigEnumInput> const &adaptable) const {
  //     return adaptable.configuration;
  //   }
  //
  //   Supercell const &
  //   Adapter<Supercell, ConfigEnumData<ConfigEnumStrain, ConfigEnumInput>>::operator()(
  //     ConfigEnumData<ConfigEnumStrain, ConfigEnumInput> const &adaptable) const {
  //     return adaptable.configuration.supercell();
  //   }
  //
  // }

  namespace ConfigEnumIO {
    /// Template specialization to get current normal coordinate from enumerator.
    template<>
    Eigen::VectorXd get_normal_coordinate(ConfigEnumStrain const &enumerator) {
      return enumerator.normal_coordinate();
    }

    template<typename ConfigEnumDataType>
    GenericDatumFormatter<Index, ConfigEnumDataType> subwedge_index() {
      return GenericDatumFormatter<Index, ConfigEnumDataType>(
               "subwedge_index",
               "Subwedge index (determines axes defining normal coordinate).",
      [](ConfigEnumDataType const & data) {
        return data.enumerator.subwedge_index();
      });
    }
  }

  std::string ConfigEnumStrainInterface::desc() const {

    std::string description =
      "Generate strain perturbations of one or more initial configurations.\n\n";

    std::string custom_options =
      "  axes: matrix or JSON object (optional, default=identity matrix of DoF space dimension) \n\n"
      "    Coordinate axes of the DoF grid. Each element in an axis vector correponds to an \n"
      "    individual DoF. Each axis vector corresponds to a normal mode. Use the option \n"
      "    `\"print_dof_space_and_quit\": true` to print DoF space information with a glossary \n"
      "    describing which DoF is specified by which vector element. The 'axes' may be rank \n"
      "    deficient indicating enumeration should occur in a subspace of the full DoF space \n"
      "    specified by the \"dof\" value and initial enumeration state.\n\n"

      "    Example if matrix (row vector matix): \n"
      "      \"axes\" : [ \n"
      "        [1, 1, 1, 1, 1, 1], \n"
      "        [1,-1, 0,-1, 1, 0], \n"
      "        [1,-1, 0, 1,-1, 0]  \n"
      "      ] \n\n"

      "    Example if JSON object (named axis vectors): \n"
      "        \"axes\": { \n"
      "          \"q1\": [1, 1, 1, 1, 1, 1], \n"
      "          \"q2\": [1,-1, 0,-1, 1, 0], \n"
      "          \"q3\": [1,-1, 0, 1,-1, 0]  \n"
      "        } \n\n"

      "    Note: \n"
      "    - If some \"qi\" in the range [1, DoF space dimension] are missing, then enumeration \n"
      "      is performed in the subspace specified by the axes that are provided. \n\n"

      "  sym_axes: bool (optional, default=false)\n"
      "    If true, overrides \"axes\" field and instead constructs symmetry-adapted grid axes\n"
      "    as the symmetry-adapted DoF order parameters of 'config'. Run with option \n"
      "    `\"print_dof_space_and_quit\": true` to obtain the analysis report including the \n"
      "    symmetry-adapted axes.\n\n"

      "  print_dof_space_and_quit: boolean (optional, default=false) \n"
      "    If true, print DoF space information for each initial enumeration state and quit. If \n"
      "    `\"sym_axes\": true`, will also print irreducible subspaces and symmetry-adapted axes. \n\n"

      "  min: number, or array of numbers (optional, default = [0,...,0]) \n"
      "    Minimum, starting value of grid counter\n"
      "    If number, specifies using a constant array of DoF space dimension with that given value.\n"
      "    Ex: \"min\" : -0.1  ( -->  [-0.1, -0.1, ..., -0.1])\n"
      "    If array, dimension must be equal to the \"axes\" dimension.\n"
      "    Ex: \"min\" : [-0.05, -0.1, -0.1]\n\n"

      "  max: number, or array of numbers (required) \n"
      "    Maximum, final value of grid counter\n"
      "    If number, specifies using a constant array of DoF space dimension with that given value.\n"
      "    Ex: \"max\" : 0.1  ( -->  [0.1, 0.1, ..., 0.1])\n"
      "    If array, dimension must be equal to the \"axes\" dimension.\n"
      "    Ex: \"max\" : [0.05, 0.1, 0.1]\n\n"

      "  increment: number, or array of numbers (required) \n"
      "    Amount by which to increment counter elements\n"
      "    If number, specifies using a constant array of DoF space dimension with that given value.\n"
      "    Ex: \"increment\" : 0.01  ( -->  [0.01, 0.01, ..., 0.01])\n"
      "    If array, dimension must be equal to the \"axes\" dimension.\n"
      "    Ex: \"max\" : [0.005, 0.01, 0.01]\n\n"

      "  trim_corners: bool (optional, default=true)\n"
      "    If true, any grid points outside the largest ellipsoid inscribed within the extrema\n"
      "    of the grid will be discarded.\n\n";


    std::string examples =
      "  Examples:\n"
      "    To enumerate all strain perturbations of a particular configuration:\n"
      "      casm enum --method ConfigEnumStrain -i \n"
      "      '{ \n"
      "         \"confignames\" : [\"SCEL4_1_4_1_0_0_0/3\"],\n"
      "         \"increment\" : [0.01, 0.01, 0.01, 0., 0., 0.],\n"
      "         \"max\" : [0.05, 0.05, 0,05, 0., 0., 0.]\n"
      "       }' \n\n";

    return name() + ": \n\n" + description + custom_options + standard_ConfigEnumInput_help() + examples;
  }

  std::string ConfigEnumStrainInterface::name() const {
    return ConfigEnumStrain::enumerator_name;
  }

  namespace ConfigEnumStrainInterface_impl {

    // Functor to construct ConfigEnumStrain given an initial state
    struct MakeEnumerator {

      MakeEnumerator(
        ConfigEnumStrainParams const &_params,
        Eigen::MatrixXd const &_axes,
        bool _make_symmetry_adapted_axes):
        log(CASM::log()),
        params(_params),
        axes(_axes),
        make_symmetry_adapted_axes(_make_symmetry_adapted_axes) {}

      Log &log;
      ConfigEnumStrainParams const &params;
      Eigen::MatrixXd const &axes;
      bool make_symmetry_adapted_axes;

      ConfigEnumStrain operator()(Index index, std::string name, ConfigEnumInput const &initial_state) const {

        if(make_symmetry_adapted_axes) { // if sym_axes==true, make and use symmetry adapted axes

          log << "Performing DoF analysis" << std::endl;
          DoFSpace dof_space {initial_state, params.dof, axes};
          bool calc_wedges = true;
          std::vector<PermuteIterator> group = initial_state.configuration().factor_group();
          VectorSpaceSymReport sym_report = vector_space_sym_report(dof_space,
                                                                    group.begin(),
                                                                    group.end(),
                                                                    calc_wedges);

          std::string filename = "dof_analysis_" + params.dof + "." + std::to_string(index) + ".json";
          log << "Writing DoF analysis: " << filename << std::endl;
          fs::ofstream file {filename};
          jsonParser json;
          to_json(sym_report, json);
          to_json(dof_space, json, name);
          file << json << std::endl;

          ConfigEnumStrainParams symmetry_adapted_params = params;
          symmetry_adapted_params.wedges = sym_report.irreducible_wedge;

          return ConfigEnumStrain {initial_state, symmetry_adapted_params};
        }
        else {
          ConfigEnumStrainParams tmp_params = params;
          tmp_params.wedges.clear();
          tmp_params.wedges.push_back(SymRepTools::SubWedge::make_dummy(axes));

          return ConfigEnumStrain {initial_state, tmp_params};
        }
      }
    };
  }

  /// Parse the ConfigEnumStrainParams (except parse "axes" instead of wedges) from JSON input
  ///
  /// \param parser InputParser will contain ConfigEnumStrainParams if successful and error messages
  ///               otherwise
  /// \param axes Eigen::MatriXd speciying a custom user subspace / normal mode axes if provied,
  ///        otherwise set to Identity matrix of size dof_space_dimension x dof_space_dimension.
  /// \param dof_space_dimension DoF space dimension detected from initial enumeration states
  void parse(
    InputParser<ConfigEnumStrainParams> &parser,
    ConfigEnumInput const &initial_state,
    Eigen::MatrixXd &axes,
    bool sym_axes_option) {

    parser.value = notstd::make_unique<ConfigEnumStrainParams>();
    auto &params = *parser.value;

    // 1) get DoF type -------------------------------------------------
    // "dof" -> params.dof
    try {
      params.dof = xtal::get_strain_dof_key(initial_state.configuration().supercell().prim());
    }
    catch(std::exception &e) {
      parser.error.insert(e.what());
      return;
    }

    Index dof_space_dimension = get_dof_space_dimension(params.dof,
                                                        initial_state.configuration(),
                                                        initial_state.sites());

    // 2) get axes and normal coordinate grid ------------------------------------------
    auto grid_parser = parser.parse_as<AxesCounterParams>(dof_space_dimension);
    if(grid_parser->valid()) {
      axes = grid_parser->value->axes;
      params.min_val = grid_parser->value->min_val;
      params.max_val = grid_parser->value->max_val;
      params.inc_val = grid_parser->value->inc_val;
    }

    // 3) set auto_range option  ---------------------------------------
    if(!parser.self.contains("min") && sym_axes_option == true) {
      params.auto_range = true;
    }
    else {
      params.auto_range = false;
    }

    // 4) get trim_corners option  -------------------------------------
    // "trim_corners" -> params.trim_corners (bool, default true)
    parser.optional_else(params.trim_corners, "trim_corners", true);

  }

  void ConfigEnumStrainInterface::run(
    PrimClex &primclex,
    jsonParser const &json_options,
    jsonParser const &cli_options_as_json) const {

    Log &log = CASM::log();

    log.subsection().begin("ConfigEnumStrain");
    log.indent() << "Input from JSON (--input or --setings):\n" << json_options << std::endl << std::endl;
    log.indent() << "Input from `casm enum` options:\n" << cli_options_as_json << std::endl << std::endl;

    // combine JSON options and CLI options
    jsonParser json_combined = combine_configuration_enum_json_options(
                                 json_options,
                                 cli_options_as_json);

    log.indent() << "Combined Input:\n" << json_combined << std::endl << std::endl;

    // Read input data from JSON
    ParentInputParser parser {json_combined};
    std::runtime_error error_if_invalid {"Error reading ConfigEnumStrain JSON input"};

    log.custom("Checking input");

    // 1) Parse initial enumeration states ------------------
    typedef std::vector<std::pair<std::string, ConfigEnumInput>> NamedInitialEnumerationStates;
    auto input_parser_ptr = parser.parse_as<NamedInitialEnumerationStates>(
                              primclex.shared_prim(),
                              &primclex,
                              primclex.db<Supercell>(),
                              primclex.db<Configuration>());
    report_and_throw_if_invalid(parser, log, error_if_invalid);
    auto const &named_initial_states = *input_parser_ptr->value;
    log.indent() << "# of initial enumeration states: " << named_initial_states.size() << std::endl;

    // 1b) If "verbose", print
    log.subsection().begin_section<Log::verbose>();
    log.indent() << "initial enumeration states:" << std::endl;
    log.increase_indent();
    for(auto const &named_initial_state : named_initial_states) {
      log.indent() << named_initial_state.first << std::endl;
    }
    log.decrease_indent();
    log.end_section();

    // 2) Parse ConfigEnumStrainParams ------------------

    // 2a) check for "sym_axes" option:
    bool sym_axes_option;
    parser.optional_else(sym_axes_option, "sym_axes", false);
    log.indent() << "sym_axes: " << sym_axes_option << std::endl;

    // 2b) parse ConfigEnumStrainParams (except parse "axes" instead of "wedges")

    // axes: column vector matrix (user input, else Identity matrix of dimension == strain space dimension
    // - If sym_axes==false: normal mode coordinates
    // - If sym_axes==true: subspace to be symmetrized
    Eigen::MatrixXd axes;
    auto params_parser_ptr = parser.parse_as<ConfigEnumStrainParams>(named_initial_states[0].second,
                                                                     axes,
                                                                     sym_axes_option);
    report_and_throw_if_invalid(parser, log, error_if_invalid);
    ConfigEnumStrainParams const &params = *params_parser_ptr->value;
    log.indent() << "axes: \n" << axes << std::endl;
    log.indent() << "min: " << params.min_val.transpose() << std::endl;
    log.indent() << "max: " << params.max_val.transpose() << std::endl;
    log.indent() << "increment: " << params.inc_val.transpose() << std::endl;
    log.indent() << "auto_range: " << params.auto_range << std::endl;
    log.indent() << "trim_corners: " << params.trim_corners << std::endl;

    // 2c) check for "print_dof_space_and_quit" option:
    bool print_dof_space_and_quit_option;
    parser.optional_else(print_dof_space_and_quit_option, "print_dof_space_and_quit", false);
    log.indent() << "print_dof_space_and_quit: " << print_dof_space_and_quit_option << std::endl;

    if(print_dof_space_and_quit_option) {
      log.begin<Log::debug>("Print DoF Space and Quit Option");
      for(auto const &named_initial_state : named_initial_states) {
        auto const &name = named_initial_state.first;
        auto const &initial_state = named_initial_state.second;
        DoFSpace dof_space {initial_state, params.dof, axes};
        bool calc_wedges = true;
        std::vector<PermuteIterator> group = initial_state.configuration().factor_group();
        print_dof_space(log, name, dof_space, group.begin(), group.end(), sym_axes_option, calc_wedges);
      }

      log.end_section();
      return;
    }

    // 3) Parse ConfigEnumOptions ------------------
    auto options_parser_ptr = parser.parse_as<ConfigEnumOptions>(
                                ConfigEnumStrain::enumerator_name,
                                primclex,
                                primclex.settings().query_handler<Configuration>().dict());
    report_and_throw_if_invalid(parser, log, error_if_invalid);
    ConfigEnumOptions const &options = *options_parser_ptr->value;
    log.indent() << "pritive_only: " << options.primitive_only << std::endl;
    log.indent() << "filter: " << static_cast<bool>(options.filter) << std::endl;
    if(options.filter) {
      std::string filter_expression;
      parser.self.get_if(filter_expression, "filter");
      log.indent() << "filter expression: " << filter_expression << std::endl;
    }
    log.indent() << "verbosity: " << options.verbosity << std::endl;
    log.indent() << "dry_run: " << options.dry_run << std::endl;

    // 4) Enumerate configurations ------------------
    log << std::endl;
    log.begin("ConfigEnumStrain enumeration");

    ConfigEnumStrainInterface_impl::MakeEnumerator make_enumerator_f {params, axes, sym_axes_option};

    typedef ConfigEnumData<ConfigEnumStrain, ConfigEnumInput> ConfigEnumDataType;

    int sep = 2;
    int width = 12;
    std::string comment = "#";
    DataFormatter<ConfigEnumDataType> formatter {sep, width, comment};

    std::string prim_strain_metric = xtal::get_strain_metric(params.dof);
    formatter.push_back(
      ConfigEnumIO::name<ConfigEnumDataType>(),
      ConfigEnumIO::selected<ConfigEnumDataType>(),
      ConfigEnumIO::is_new<ConfigEnumDataType>(),
      ConfigEnumIO::is_existing<ConfigEnumDataType>(),
      ConfigEnumIO::is_excluded_by_filter<ConfigEnumDataType>(),
      ConfigEnumIO::initial_state_index<ConfigEnumDataType>(),
      ConfigEnumIO::initial_state_name<ConfigEnumDataType>(),
      ConfigEnumIO::initial_state_configname<ConfigEnumDataType>(),
      ConfigEnumIO::subwedge_index<ConfigEnumDataType>(),
      ConfigEnumIO::normal_coordinate<ConfigEnumDataType>(),
      make_datum_formatter_adapter<ConfigEnumDataType, Configuration>(
        ConfigIO::DoFStrain(prim_strain_metric))
    );
    if(prim_strain_metric != "F") { // is this ever not the case?
      formatter.push_back(
        make_datum_formatter_adapter<ConfigEnumDataType, Configuration>(ConfigIO::DoFStrain("F"))
      );
    }

    enumerate_configurations(
      primclex,
      options,
      make_enumerator_f,
      named_initial_states.begin(),
      named_initial_states.end(),
      formatter);

    log.indent() << "enumeration complete" << std::endl << std::endl;
    log.end_section();
  }

}
