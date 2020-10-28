#include "casm/app/enum/enumerate_configurations_impl.hh"
#include "casm/app/enum/methods/ConfigEnumSiteDoFsInterface.hh"
#include "casm/clex/ConfigEnumSiteDoFs.hh"

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

  std::string ConfigEnumSiteDoFsInterface::desc() const {

    std::string description =

      "The ConfigEnumSiteDoFs method generates configurations corresponding excitations of collective \n"
      "site DoF (normal) modes. The input consists of three components:\n\n"

      "1) Specification of one or more site DoF spaces to sample. This is done via:\n\n"

      "   - Choice of site DoF type\n"
      "   - Specification of one or more configurations.\n"
      "   - Selection of particular sites in the configuration(s). Sites may be selected by indicating\n"
      "     particular sites (\"sites\"), cluster of sites (\"cluster_specs\"), or by all sites on\n"
      "     particular sublattices (\"sublats\"). By default, if none are specified, then all sites in\n"
      "     the configuration are selected.\n\n"

      "   The total dimension of the specified site DoF space is equal to the sum over all chosen\n"
      "   sites of the site DoF axes dimension. If the site DoF axes are the same on all selected\n"
      "   sites, then the total dimension is (number of sites x site DoF dimension). If the site\n"
      "   DoF space is restricted on some sublattices, then the total dimension may be smaller.\n\n"

      "     Examples:\n"
      "     - The standard \"disp\" basis is used on all sublattices. If 4 sites are selected:\n"
      "       - The total site DoF space dimension is 12 = (4 sites) * (3 dimensional disp)\n"
      "     - The standard \"disp\" basis is used on sublattice b=0, and on sublattice b=1 the\n"
      "       \"disp\" axes are set to only allow 1 dimensional displacements. If 2 sites from\n"
      "       sublattice b=0 are selected and 1 site from the sublattice b=1 is selected:\n"
      "       - The total site DoF space dimension is 7 =\n"
      "           (2 sites from sublattice b=0) * (3 dimensional disp) +\n"
      "           (1 site from sublatice b=1) * (1 dimensional disp)\n\n"

      "     Notes on \"cluster_specs\":\n"
      "     - The \"cluster_specs\" option may not be used with the \"sites\" or \"sublats\" options.\n"
      "     - The cluster orbits are generated using the configuration factor group symmetry and\n"
      "       then any orbits that are duplicated under periodic boundary conditions are removed.\n\n\n"


      "2) Specification of a normal modes in the site DoF space. This can be done through a\n"
      "   combination of the \"axes\" and \"sym_axes\" options:\n\n"

      "   - The parameter \"axes\" is an optional row matrix of normal coordinate axes. If it is not\n"
      "     provided, then it is set to the identity matrix with dimension equal to the dimension of\n"
      "     the site DoF space specified by the choice of configuration and selected sites. It is not\n"
      "     required to be full rank (i.e. number of axes rows < number of axes columns is valid), in\n"
      "     which case it means Configurations are generated in a subspace.\n\n"

      "     Notes on \"axes\":\n"
      "     - Each column in \"axes\" corresponds to an individual prim DoF, which is printed to screen\n"
      "       - Examples:\n"
      "         - Column c: type=\"disp\", component=(\"dx\", [1, 0, 0]), sublattice index=b, unit cell=(i,j,k)\n"
      "         - Column c: type=\"disp\", component=(\"dxy\", [1, 1, 0]), sublattice index=b, unit cell=(i,j,k)\n"
      "     - Each row is an axis in the total site DoF space\n\n"

      "   - If the optional parameter \"sym_axes\" is true (default=false), then CASM will generate\n"
      "     symmetry adapted normal coordinate axes in the space (may be a subspace) specified by\n"
      "     \"axes\" for each initial enumeration state. This means that CASM:\n"
      "     - Finds the configuration factor group (symmetry operations that keep the supercell\n"
      "       lattice invariant and the configuration DoF values invariant)\n"
      "     - Finds the subgroup which also keeps the selected sites invariant (does not permute\n"
      "       selected sites with unselected sites)\n"
      "     - Calculates the irreducible subspaces of the site DoF space under that subgroup\n"
      "     - Uses the axes of the irreducible subspaces as the normal coordinate axes\n\n"

      "     Notes on \"sym_axes\":\n"
      "     - If \"sym_axes\"==false, the coordinate axes are used directly for the normal modes\n"
      "     - If \"sym_axes\"==true, symmetry adapted normal modes are generated in the subspace defined\n"
      "       by \"axes\" (default is total site DoF space) and printed to screen\n"
      "     - The symmetry adapted axes can also be calculated via the `casm sym` command\n"
      "     - The user may take the symmetry adapted axes, rotate the irreducible subspaces, and use\n"
      "       that as the \"axes\" input, with \"sym_axes\"=false, to customize the choice of normal\n"
      "       coordinates.\n\n\n"


      "3) Choice of linear combinations of normal modes to apply to the chosen Configuration:\n\n"

      "   Even if \"axes\" are rank deficient, the site DoF space defined by axes may quickly become\n"
      "   very high dimensional (number of sites x mean site DoF dimension), so rather than sample the\n"
      "   entire space, ConfigEnumSiteDoFs perturbs the input configuration by applying a linear\n"
      "   combination of normal modes.\n\n"

      "   The amplitudes of the normal modes is specified with the \"min\", \"increment\", and \"max\"\n"
      "   parameters. These may be scalar valued, to set sampled amplitudes to be the same along each\n"
      "   normal coordinate axes dimension. Of they may be vector valued in order to customize the\n"
      "   sampled amplitudes along different dimensions. If the total dimension of the site DoF\n"
      "   varies with choice of input configurations and selected sites, then only the scalar input\n"
      "   option is allowed.\n\n"

      "   The parameters \"min_nonzero\" and \"max_nonzero\" specifies how many normal mode amplitudes\n"
      "   should be nonzero (inclusive range [min_nonzero, max_nonzero]). The method generates\n"
      "   all n choose k (n=site DoF space dimension, k iterates through [min_nonzer, max_nonzero])\n"
      "   combinations of normal modes in that range, and for each combination applies all the\n"
      "   k chosen normal modes with amplitudes specified by \"min\" / \"increment\" / \"max\". Note that\n"
      "   this may quickly become very large, depending on n, k, and the range specified by \"min\" /\n"
      "   \"increment\" / \"max\".\n\n";

    std::string custom_options =

      "  dof: string (required) \n"
      "    Name of degree of freedom for which normal coordinates are to be generated.\n"
      "    Must be one of the degrees of freedom under consideration in the current project,\n"
      "    as specified in prim.json.\n\n"

      "  axes: matrix or JSON object (optional, default=identity matrix of DoF space dimension) \n\n"
      "    Coordinate axes of the DoF grid. Each element in an axis vector correponds to an \n"
      "    individual DoF. Each axis vector corresponds to a normal mode. Use the option \n"
      "    `\"print_dof_space_and_quit\": true` to print DoF space information with a glossary \n"
      "    describing which DoF is specified by which vector element. The 'axes' may be rank \n"
      "    deficient indicating enumeration should occur in a subspace of the full DoF space \n"
      "    specified by the \"dof\" value and initial enumeration state.\n"

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

      "      Note: \n"
      "      - If some \"qi\" in the range [1, DoF space dimension] are missing, then enumeration \n"
      "        is performed in the subspace specified by the axes that are provided. \n\n"

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

      "  min_nonzero: integer (optional, default = 0) \n"
      "    Minimum number of coordinate amplitudes that are allowed\n"
      "    to be nonzero. Must be less than or equal to the \"axes\" dimension.\n\n"

      "  max_nonzero: integer (optional, default = axes.rows()) \n"
      "    Maximum number of coordinate amplitudes that are allowed\n"
      "    to be nonzero. Must be less than or equal to the \"axes\" dimension.\n\n";

    std::string examples =
      "  Examples:\n"
      "    To enumerate all DoF perturbations of a particular configuration:\n"
      "      casm enum --method ConfigEnumSiteDoFs -i \n"
      "      '{ \n"
      "        \"confignames\": [ \n"
      "          \"SCEL4_1_4_1_0_0_0/3\"\n"
      "        ],\n"
      "        \"max\": 0.21, \n"
      "        \"increment\": 0.05, \n"
      "        } \n"
      "      }' \n\n";

    return name() + ": \n\n" + description + custom_options + standard_ConfigEnumInput_help() + examples;
  }

  std::string ConfigEnumSiteDoFsInterface::name() const {
    return ConfigEnumSiteDoFs::enumerator_name;
  }

  /// Parse the "min" / "increment" / "max" scalar or vector values
  void parse_counter_value(
    InputParser<ConfigEnumSiteDoFsParams> &parser,
    Eigen::VectorXd &counter_value,
    std::string attribute_name,
    Index n_axes_columns,
    bool has_default_value,
    Eigen::VectorXd const &default_value) {

    if(parser.self.contains(attribute_name) && parser.self[attribute_name].is_number()) {
      counter_value = Eigen::VectorXd::Constant(
                        n_axes_columns,
                        parser.self[attribute_name].get<double>());
    }
    else if(parser.self.contains(attribute_name) && parser.self[attribute_name].is_array()) {
      parser.optional(counter_value, attribute_name);
      if(counter_value.size() != n_axes_columns) {
        std::stringstream msg;
        // note: message refers to JSON axes rows, but that means params.axes.cols()
        msg << "Error: \"" << attribute_name << "\" array size must match dimension of \"axes\" rows";
        parser.error.insert(msg.str());
      }
    }
    else if(has_default_value) {
      counter_value = default_value;
    }
    else {
      std::stringstream msg;
      msg << "Error: \"" << attribute_name << "\" must be a number or array of number";
      parser.error.insert(msg.str());
    }
  }

  /// Parse the ConfigEnumSiteDoFsParams JSON input
  void parse(
    InputParser<ConfigEnumSiteDoFsParams> &parser,
    ConfigEnumInput const &initial_state) {

    parser.value = notstd::make_unique<ConfigEnumSiteDoFsParams>();
    auto &params = *parser.value;

    // 1) get DoF type -------------------------------------------------
    // "dof" -> params.dof
    parser.require(params.dof, "dof");

    Index dof_space_dimension = get_dof_space_dimension(params.dof,
                                                        initial_state.configuration(),
                                                        initial_state.sites());

    // 2) get axes -----------------------------------------------------
    // "axes".transpose() -> params.axes  (transpose of each other)
    Eigen::MatrixXd default_axes = Eigen::MatrixXd::Identity(dof_space_dimension, dof_space_dimension);
    Eigen::MatrixXd row_vector_axes;
    parser.optional_else(row_vector_axes, "axes", default_axes);
    params.axes = row_vector_axes.transpose();

    // check axes dimensions:
    if(params.axes.rows() != dof_space_dimension) {
      // Note: message "columns" refers to JSON input axes, transpose of params.axes
      std::stringstream msg;
      msg << "Number of columns of \"axes\" must be equal to site DoF space dimension ("
          << dof_space_dimension << "). Size as parsed: " << params.axes.rows();
      parser.error.insert(msg.str());
    }
    if(params.axes.cols() > dof_space_dimension) {
      // Note: message "rows" refers to JSON input axes, transpose of params.axes
      std::stringstream msg;
      msg << "Number of coordinate axes (number of rows of \"axes\") must be less than or equal to "
          "site DoF space dimension (" << dof_space_dimension << "). Number of axes parsed: "
          << params.axes.cols();
      parser.error.insert(msg.str());
    }


    // 3) get normal coordinate grid -----------------------------------

    // "min" -> params.min_val  (number or array of number, else zeros vector)
    Eigen::VectorXd default_value = Eigen::VectorXd::Zero(params.axes.cols());
    params.min_val = parse_vector_from_number_or_array(parser, "min", params.axes.cols(), &default_value);
    // note that help indicates default==axes.rows(), but that is params.axes.cols()

    // "max" -> params.max_val (number or array of number, required)
    params.max_val = parse_vector_from_number_or_array(parser, "max", params.axes.cols());

    // "increment" -> params.inc_val (number or array of number, required)
    params.inc_val = parse_vector_from_number_or_array(parser, "increment", params.axes.cols());


    // 4) get min/max nonzero amplitudes -----------------------------------

    // "min_nonzero" -> params.min_nonzero
    parser.optional_else(params.min_nonzero, "min_nonzero", Index {0});

    // "max_nonzero" -> params.max_nonzero
    // note that help indicates default==axes.rows(), but that is params.axes.cols()
    parser.optional_else(params.min_nonzero, "max_nonzero", Index {params.axes.cols()});
  }

  void require_all_input_have_the_same_number_of_selected_sites(
    InputParser<std::vector<std::pair<std::string, ConfigEnumInput>>> &parser) {

    // TODO: should this check if the supercells & sites are the same, and not just the same number?
    // if sym_axes==false, then using different sites / supercells might be bad

    if(parser.valid() && parser.value != nullptr) {
      auto const &named_initial_states = *parser.value;
      Index nsites = named_initial_states[0].second.sites().size();
      for(auto const &name_value_pair : named_initial_states) {
        if(name_value_pair.second.sites().size() != nsites) {
          std::stringstream msg;
          msg << "Error in ConfigEnumSiteDoFs: Starting configurations or supercells passed to "
              "must all have the same number of selected sites.";
          parser.error.insert(msg.str());
        }
      }
    }
  }

  namespace ConfigEnumSiteDoFsInterface_impl {

    // This functor constructs a ConfigEnumSiteDoFs enumerator for a given initial state
    struct MakeEnumerator {

      MakeEnumerator(
        ConfigEnumSiteDoFsParams const &_params,
        bool _make_symmetry_adapted_axes):
        log(CASM::log()),
        params(_params),
        make_symmetry_adapted_axes(_make_symmetry_adapted_axes) {}

      Log &log;
      ConfigEnumSiteDoFsParams const &params;
      bool make_symmetry_adapted_axes;

      ConfigEnumSiteDoFs operator()(std::string name, ConfigEnumInput const &initial_state) const {

        if(make_symmetry_adapted_axes) { // if sym_axes==true, make and use symmetry adapted axes

          log.begin(std::string("DoF Vector Space Symmetry Report: ") + name);
          log << "For large spaces this may be slow..." << std::endl;

          DoFSpace dof_space {initial_state, params.dof, params.axes};
          bool calc_wedges = false;
          std::vector<PermuteIterator> group = initial_state.configuration().factor_group();
          VectorSpaceSymReport sym_report = vector_space_sym_report(dof_space,
                                                                    group.begin(),
                                                                    group.end(),
                                                                    calc_wedges);
          log << jsonParser {sym_report} << std::endl;
          log.end_section();

          ConfigEnumSiteDoFsParams symmetry_adapted_params = params;
          symmetry_adapted_params.axes = sym_report.symmetry_adapted_dof_subspace;

          return ConfigEnumSiteDoFs {initial_state, symmetry_adapted_params};
        }
        else { // if sym_axes==false, use input or default axes

          return ConfigEnumSiteDoFs {initial_state, params};
        }
      }
    };
  }

  void ConfigEnumSiteDoFsInterface::run(
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
    std::runtime_error error_if_invalid {"Error reading ConfigEnumSiteDoFs JSON input"};

    // 1) Parse initial enumeration states ------------------
    typedef std::vector<std::pair<std::string, ConfigEnumInput>> NamedInitialEnumerationStates;
    auto input_parser_ptr = parser.parse_as<NamedInitialEnumerationStates>(
                              primclex.shared_prim(),
                              &primclex,
                              primclex.db<Supercell>(),
                              primclex.db<Configuration>());
    require_all_input_have_the_same_number_of_selected_sites(*input_parser_ptr);
    report_and_throw_if_invalid(parser, log, error_if_invalid);
    auto const &named_initial_states = *input_parser_ptr->value;

    // 2) Parse ConfigEnumSiteDoFsParams ------------------
    auto params_parser_ptr = parser.parse_as<ConfigEnumSiteDoFsParams>(named_initial_states[0].second);
    report_and_throw_if_invalid(parser, log, error_if_invalid);
    ConfigEnumSiteDoFsParams const &params = *params_parser_ptr->value;

    // 2b) check for "sym_axes" option:
    bool sym_axes_option;
    parser.optional_else(sym_axes_option, "sym_axes", false);

    // 2c) check for "print_dof_space_and_quit" option:
    bool print_dof_space_and_quit_option;
    parser.optional_else(print_dof_space_and_quit_option, "print_dof_space_and_quit", false);

    if(print_dof_space_and_quit_option) {
      for(auto const &named_initial_state : named_initial_states) {
        auto const &name = named_initial_state.first;
        auto const &initial_state = named_initial_state.second;
        DoFSpace dof_space {initial_state, params.dof, params.axes};
        bool calc_wedges = false;
        std::vector<PermuteIterator> group = initial_state.configuration().factor_group();
        print_dof_space(log, name, dof_space, group.begin(), group.end(), sym_axes_option, calc_wedges);
      }
      return;
    }

    // 3) Parse EnumerateConfigurationsOptions ------------------
    auto options_parser_ptr = parser.parse_as<EnumerateConfigurationsOptions>(
                                ConfigEnumSiteDoFs::enumerator_name,
                                primclex,
                                primclex.settings().query_handler<Configuration>().dict());
    report_and_throw_if_invalid(parser, log, error_if_invalid);
    EnumerateConfigurationsOptions const &options = *options_parser_ptr->value;

    // 4) Enumerate configurations ------------------

    ConfigEnumSiteDoFsInterface_impl::MakeEnumerator make_enumerator_f {params, sym_axes_option};

    enumerate_configurations(
      options,
      make_enumerator_f,
      named_initial_states.begin(),
      named_initial_states.end(),
      primclex.db<Supercell>(),
      primclex.db<Configuration>());
  }

}
