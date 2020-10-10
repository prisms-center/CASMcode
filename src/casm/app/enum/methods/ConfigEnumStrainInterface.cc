#include "casm/app/enum/enumerate_configurations_impl.hh"
#include "casm/app/enum/methods/ConfigEnumStrainInterface.hh"
#include "casm/clex/ConfigEnumStrain.hh"

#include "casm/app/APICommand.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler_impl.hh"
#include "casm/app/enum/standard_ConfigEnumInput_help.hh"
#include "casm/app/enum/io/enumerate_configurations_json_io.hh"
#include "casm/app/enum/io/json_io.hh"
#include "casm/app/enum/io/stream_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/DoFSpace.hh"
#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"
#include "casm/symmetry/SymRepTools.hh"
#include "casm/symmetry/io/json/SymRepTools.hh"

namespace CASM {

  std::string ConfigEnumStrainInterface::desc() const {

    std::string description =
      "Generate strain perturbations of one or more initial configurations.\n\n";

    std::string custom_options =
      "  axes: matrix of doubles (optional, default=identity matrix of DoF space dimension) \n"
      "    Coordinate axes of the DoF grid. Each column correponds to an individual DoF. Each row \n"
      "    corresponds to a normal mode. Use the option `\"print_dof_space_and_quit\": true` to print \n"
      "    DoF space information with a glossary describing which DoF is specified by which column \n"
      "    for a particular initial enumeration state. The 'axes' matrix may be rank deficient \n"
      "    indicating enumeration should occur in a subspace of the full DoF space specified by the \n"
      "    \"dof\" value and initial enumeration state.\n"
      "    Ex: \"axes\" : [[1, 1, 1, 1, 1, 1],\n"
      "                  [1,-1, 0,-1, 1, 0],\n"
      "                  [1,-1, 0, 1,-1, 0]]\n\n"

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
      "    If array, dimension must be equal to DoF space dimension.\n"
      "    Ex: \"min\" : [-0.05, -0.1, -0.1]\n\n"

      "  max: number, or array of numbers (required) \n"
      "    Maximum, final value of grid counter\n"
      "    If number, specifies using a constant array of DoF space dimension with that given value.\n"
      "    Ex: \"max\" : 0.1  ( -->  [0.1, 0.1, ..., 0.1])\n"
      "    If array, dimension must be equal to DoF space dimension.\n"
      "    Ex: \"max\" : [0.05, 0.1, 0.1]\n\n"

      "  increment: number, or array of numbers (required) \n"
      "    Amount by which to increment counter elements\n"
      "    If number, specifies using a constant array of DoF space dimension with that given value.\n"
      "    Ex: \"increment\" : 0.01  ( -->  [0.01, 0.01, ..., 0.01])\n"
      "    If array, dimension must be equal to DoF space dimension.\n"
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

      ConfigEnumStrain operator()(std::string name, ConfigEnumInput const &initial_state) const {

        if(make_symmetry_adapted_axes) { // if sym_axes==true, make and use symmetry adapted axes

          log.begin(std::string("DoF Vector Space Symmetry Report: ") + name);

          DoFSpace dof_space {initial_state, params.dof, axes};
          bool calc_wedges = true;
          std::vector<PermuteIterator> group = initial_state.configuration().factor_group();
          VectorSpaceSymReport sym_report = vector_space_sym_report(dof_space,
                                                                    group.begin(),
                                                                    group.end(),
                                                                    calc_wedges);
          log << jsonParser {sym_report} << std::endl;
          log.end_section();

          ConfigEnumStrainParams symmetry_adapted_params = params;
          symmetry_adapted_params.wedges = sym_report.irreducible_wedge;

          return ConfigEnumStrain {initial_state, symmetry_adapted_params};
        }
        else {
          ConfigEnumStrainParams tmp_params = params;
          tmp_params.wedges.clear();
          tmp_params.wedges.push_back(SymRepTools::SubWedge::make_dummy(axes));

          return ConfigEnumStrain {initial_state, params};
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
    parser.require(params.dof, "dof");

    Index dof_space_dimension = get_dof_space_dimension(params.dof,
                                                        initial_state.configuration(),
                                                        initial_state.sites());

    // 2) get axes -----------------------------------------------------
    // "axes".transpose() -> axes  (transpose of each other)
    Eigen::MatrixXd default_axes = Eigen::MatrixXd::Identity(
                                     dof_space_dimension,
                                     dof_space_dimension);
    Eigen::MatrixXd row_vector_axes;
    parser.optional_else(row_vector_axes, "axes", default_axes);
    axes = row_vector_axes.transpose();

    // check axes dimensions:
    if(axes.rows() != dof_space_dimension) {
      // Note: message "columns" refers to JSON input axes, transpose of axes
      std::stringstream msg;
      msg << "Number of columns of \"axes\" must be equal to DoF space dimension ("
          << dof_space_dimension << "). Size as parsed: " << axes.rows();
      parser.error.insert(msg.str());
    }
    if(axes.cols() > dof_space_dimension) {
      // Note: message "rows" refers to JSON input axes, transpose of axes
      std::stringstream msg;
      msg << "Number of coordinate axes (number of rows of \"axes\") must be less than or equal to "
          "DoF space dimension (" << dof_space_dimension << "). Number of axes parsed: "
          << axes.cols();
      parser.error.insert(msg.str());
    }

    // 3) get normal coordinate grid -----------------------------------

    // "min" -> params.min_val  (number or array of number, else zeros vector)
    Eigen::VectorXd default_value = Eigen::VectorXd::Zero(axes.cols());
    params.min_val = parse_vector_from_number_or_array(parser, "min", axes.cols(), &default_value);

    // "max" -> params.max_val (number or array of number, required)
    params.max_val = parse_vector_from_number_or_array(parser, "max", axes.cols());

    // "increment" -> params.inc_val (number or array of number, required)
    params.inc_val = parse_vector_from_number_or_array(parser, "increment", axes.cols());

    // 4) set auto_range option  ---------------------------------------
    if(!parser.self.contains("min") && sym_axes_option == true) {
      params.auto_range = true;
    }
    else {
      params.auto_range = false;
    }

    // 5) get trim_corners option  -------------------------------------
    // "trim_corners" -> params.trim_corners (bool, default true)
    parser.optional_else(params.trim_corners, "trim_corners", true);

  }

  void ConfigEnumStrainInterface::run(
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
    std::runtime_error error_if_invalid {"Error reading ConfigEnumStrain JSON input"};

    // 1) Parse initial enumeration states ------------------
    typedef std::vector<std::pair<std::string, ConfigEnumInput>> NamedInitialEnumerationStates;
    auto input_parser_ptr = parser.parse_as<NamedInitialEnumerationStates>(
                              primclex.shared_prim(),
                              &primclex,
                              primclex.db<Supercell>(),
                              primclex.db<Configuration>());
    auto const &named_initial_states = *input_parser_ptr->value;

    // 2) Parse ConfigEnumStrainParams ------------------

    // 2a) check for "sym_axes" option:
    bool sym_axes_option;
    parser.optional_else(sym_axes_option, "sym_axes", false);

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

    // 2c) check for "print_dof_space_and_quit" option:
    bool print_dof_space_and_quit_option;
    parser.optional_else(print_dof_space_and_quit_option, "print_dof_space_and_quit", false);

    if(print_dof_space_and_quit_option) {
      for(auto const &named_initial_state : named_initial_states) {
        auto const &name = named_initial_state.first;
        auto const &initial_state = named_initial_state.second;
        DoFSpace dof_space {initial_state, params.dof, axes};
        bool calc_wedges = true;
        std::vector<PermuteIterator> group = initial_state.configuration().factor_group();
        print_dof_space(log, name, dof_space, group.begin(), group.end(), sym_axes_option, calc_wedges);
      }
      return;
    }

    // 3) Parse EnumerateConfigurationsOptions ------------------
    auto options_parser_ptr = parser.parse_as<EnumerateConfigurationsOptions>(
                                ConfigEnumStrain::enumerator_name,
                                primclex,
                                primclex.settings().query_handler<Configuration>().dict());
    report_and_throw_if_invalid(parser, log, error_if_invalid);
    EnumerateConfigurationsOptions const &options = *options_parser_ptr->value;

    // 4) Enumerate configurations ------------------

    ConfigEnumStrainInterface_impl::MakeEnumerator make_enumerator_f {params, axes, sym_axes_option};

    enumerate_configurations(
      options,
      make_enumerator_f,
      named_initial_states.begin(),
      named_initial_states.end(),
      primclex.db<Supercell>(),
      primclex.db<Configuration>());
  }

}
