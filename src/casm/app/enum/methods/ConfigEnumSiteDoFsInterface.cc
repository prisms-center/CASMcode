#include "casm/app/enum/methods/ConfigEnumSiteDoFsInterface.hh"

#include "casm/app/APICommand.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler_impl.hh"
#include "casm/app/enum/dataformatter/ConfigEnumIO_impl.hh"
#include "casm/app/enum/enumerate_configurations_impl.hh"
#include "casm/app/enum/io/enumerate_configurations_json_io.hh"
#include "casm/app/enum/io/stream_io_impl.hh"
#include "casm/app/enum/standard_ConfigEnumInput_help.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/ConfigEnumSiteDoFs.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/DoFSpace.hh"
#include "casm/enumerator/io/dof_space_analysis.hh"
#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"
#include "casm/symmetry/SupercellSymInfo_impl.hh"
#include "casm/symmetry/SymRepTools.hh"
#include "casm/symmetry/io/json/SymRepTools.hh"

namespace CASM {

std::string ConfigEnumSiteDoFsInterface::desc() const {
  std::string description =

      "The ConfigEnumSiteDoFs method generates configurations corresponding "
      "excitations of collective \n"
      "site DoF (normal) modes. The input consists of three components:\n\n"

      "1) Specification of one or more site DoF spaces to sample. This is done "
      "via:\n\n"

      "   - Choice of site DoF type\n"
      "   - Specification of one or more configurations.\n"
      "   - Selection of particular sites in the configuration(s). Sites may "
      "be selected by indicating\n"
      "     particular sites (\"sites\"), cluster of sites "
      "(\"cluster_specs\"), or by all sites on\n"
      "     particular sublattices (\"sublats\"). By default, if none are "
      "specified, then all sites in\n"
      "     the configuration are selected.\n\n"

      "   The total dimension of the specified site DoF space is equal to the "
      "sum over all chosen\n"
      "   sites of the site DoF axes dimension. If the site DoF axes are the "
      "same on all selected\n"
      "   sites, then the total dimension is (number of sites x site DoF "
      "dimension). If the site\n"
      "   DoF space is restricted on some sublattices, then the total "
      "dimension may be smaller.\n\n"

      "     Examples:\n"
      "     - The standard \"disp\" basis is used on all sublattices. If 4 "
      "sites are selected:\n"
      "       - The total site DoF space dimension is 12 = (4 sites) * (3 "
      "dimensional disp)\n"
      "     - The standard \"disp\" basis is used on sublattice b=0, and on "
      "sublattice b=1 the\n"
      "       \"disp\" axes are set to only allow 1 dimensional displacements. "
      "If 2 sites from\n"
      "       sublattice b=0 are selected and 1 site from the sublattice b=1 "
      "is selected:\n"
      "       - The total site DoF space dimension is 7 =\n"
      "           (2 sites from sublattice b=0) * (3 dimensional disp) +\n"
      "           (1 site from sublatice b=1) * (1 dimensional disp)\n\n"

      "     Notes on \"cluster_specs\":\n"
      "     - The \"cluster_specs\" option may not be used with the \"sites\" "
      "or \"sublats\" options.\n"
      "     - The cluster orbits are generated using the configuration factor "
      "group symmetry and\n"
      "       then any orbits that are duplicated under periodic boundary "
      "conditions are removed.\n\n\n"

      "2) Specification of a normal modes in the site DoF space. This can be "
      "done through a\n"
      "   combination of the \"axes\" and \"sym_axes\" options:\n\n"

      "   - The parameter \"axes\" is an optional row matrix of normal "
      "coordinate axes. If it is not\n"
      "     provided, then it is set to the identity matrix with dimension "
      "equal to the dimension of\n"
      "     the site DoF space specified by the choice of configuration and "
      "selected sites. It is not\n"
      "     required to be full rank (i.e. number of axes rows < number of "
      "axes columns is valid), in\n"
      "     which case it means Configurations are generated in a subspace.\n\n"

      "     Notes on \"axes\":\n"
      "     - Each column in \"axes\" corresponds to an individual prim DoF, "
      "which is printed to screen\n"
      "       - Examples:\n"
      "         - Column c: type=\"disp\", component=(\"dx\", [1, 0, 0]), "
      "sublattice index=b, unit cell=(i,j,k)\n"
      "         - Column c: type=\"disp\", component=(\"dxy\", [1, 1, 0]), "
      "sublattice index=b, unit cell=(i,j,k)\n"
      "     - Each row is an axis in the total site DoF space\n\n"

      "   - If the optional parameter \"sym_axes\" is true (default=false), "
      "then CASM will generate\n"
      "     symmetry adapted normal coordinate axes in the space (may be a "
      "subspace) specified by\n"
      "     \"axes\" for each initial enumeration state. This means that "
      "CASM:\n"
      "     - Finds the configuration factor group (symmetry operations that "
      "keep the supercell\n"
      "       lattice invariant and the configuration DoF values invariant)\n"
      "     - Finds the subgroup which also keeps the selected sites invariant "
      "(does not permute\n"
      "       selected sites with unselected sites)\n"
      "     - Calculates the irreducible subspaces of the site DoF space under "
      "that subgroup\n"
      "     - Uses the axes of the irreducible subspaces as the normal "
      "coordinate axes\n\n"

      "     Notes on \"sym_axes\":\n"
      "     - If \"sym_axes\"==false, the coordinate axes are used directly "
      "for the normal modes\n"
      "     - If \"sym_axes\"==true, symmetry adapted normal modes are "
      "generated in the subspace defined\n"
      "       by \"axes\" (default is total site DoF space) and printed to "
      "screen\n"
      "     - The symmetry adapted axes can also be calculated via the `casm "
      "sym` command\n"
      "     - The user may take the symmetry adapted axes, rotate the "
      "irreducible subspaces, and use\n"
      "       that as the \"axes\" input, with \"sym_axes\"=false, to "
      "customize the choice of normal\n"
      "       coordinates.\n\n\n"

      "3) Choice of linear combinations of normal modes to apply to the chosen "
      "Configuration:\n\n"

      "   Even if \"axes\" are rank deficient, the site DoF space defined by "
      "axes may quickly become\n"
      "   very high dimensional (number of sites x mean site DoF dimension), "
      "so rather than sample the\n"
      "   entire space, ConfigEnumSiteDoFs perturbs the input configuration by "
      "applying a linear\n"
      "   combination of normal modes.\n\n"

      "   The amplitudes of the normal modes is specified with the \"min\", "
      "\"increment\", and \"max\"\n"
      "   parameters. These may be scalar valued, to set sampled amplitudes to "
      "be the same along each\n"
      "   normal coordinate axes dimension. Of they may be vector valued in "
      "order to customize the\n"
      "   sampled amplitudes along different dimensions. If the total "
      "dimension of the site DoF\n"
      "   varies with choice of input configurations and selected sites, then "
      "only the scalar input\n"
      "   option is allowed.\n\n"

      "   The parameters \"min_nonzero\" and \"max_nonzero\" specifies how "
      "many normal mode amplitudes\n"
      "   should be nonzero (inclusive range [min_nonzero, max_nonzero]). The "
      "method generates\n"
      "   all n choose k (n=site DoF space dimension, k iterates through "
      "[min_nonzer, max_nonzero])\n"
      "   combinations of normal modes in that range, and for each combination "
      "applies all the\n"
      "   k chosen normal modes with amplitudes specified by \"min\" / "
      "\"increment\" / \"max\". Note that\n"
      "   this may quickly become very large, depending on n, k, and the range "
      "specified by \"min\" /\n"
      "   \"increment\" / \"max\".\n\n";

  std::string custom_options =

      "  dof: string (required)                                           \n"
      "    Name of degree of freedom for which normal coordinates are to  \n"
      "    be generated. Must be one of the degrees of freedom under      \n"
      "    consideration in the current project, as specified in          \n"
      "    prim.json.                                                     \n\n"

      "  axes: matrix or JSON object (optional)                           \n"
      "    Coordinate axes of the DoF grid. Default value is the identity \n"
      "    matrix of DoF space dimension Each element in an axis vector   \n"
      "    correponds to an individual DoF. Each axis vector corresponds  \n"
      "    to a normal mode. Use the option                               \n"
      "    `\"print_dof_space_and_quit\": true` to print DoF space        \n"
      "    information with a glossary describing which DoF is specified  \n"
      "    by which vector element. The 'axes' may be rank deficient      \n"
      "    indicating enumeration should occur in a subspace of the full  \n"
      "    DoF space specified by the \"dof\" value and initial           \n"
      "    enumeration state.                                             \n\n"

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
      "    - If some \"qi\" in the range [1, DoF space dimension] are     \n"
      "      missing, then enumeration is performed in the subspace       \n"
      "      specified by the axes that are provided.                     \n\n"

      "  sym_axes: bool (optional, default=false)                         \n"
      "    If true, constructs symmetry-adapted grid axes as the symmetry-\n"
      "    adapted DoF order parameters of the input state in the space   \n"
      "    specified by \"axes\". Run with option                         \n"
      "    `\"print_dof_space_and_quit\": true` to obtain the analysis    \n"
      "    report including the symmetry-adapted axes before doing the    \n"
      "    enumeration. If `\"sym_axes\": true`, the enumeration grid may \n"
      "    not be specified using vector-valued \"min\", \"max\",         \n"
      "    \"increment\", or \"num\".                                     \n\n"

      "  exclude_homogeneous_modes: bool (optional, default=false)        \n"
      "   If true, homogeneous modes are removed from \"axes\" before     \n"
      "   construcing symmetry adapted axes. In the case of \"disp\" DoF, \n"
      "   this excludes all rigid translations. This option is allowed    \n"
      "   only for `\"sym_axes\": true`.                                  \n\n"

      "  print_dof_space_and_quit: boolean (optional, default=false)      \n"
      "    If true, print DoF space information for each initial          \n"
      "    enumeration state and quit. If `\"sym_axes\": true`, will also \n"
      "    print irreducible subspaces, symmetry-adapted axes, and axes   \n"
      "    excluding homogeneous modes.                                   \n\n"

      "  min: number, or array of numbers (optional)                      \n"
      "    Minimum, starting value of grid counter. If number, specifies  \n"
      "    using a constant array of DoF space dimension with that given  \n"
      "    value.                                                         \n"
      "    Ex: \"min\" : -0.1  ( -->  [-0.1, -0.1, ..., -0.1])            \n"
      "    If array, dimension must be equal to the \"axes\" dimension and\n"
      "    `\"sym_axes\" must be false`.                                  \n"
      "    Ex: \"min\" : [-0.05, -0.1, -0.1]                              \n\n"

      "  max: number, or array of numbers (required)                      \n"
      "    Maximum, final value of grid counter. If number, specifies     \n"
      "    using a constant array of DoF space dimension with that given  \n"
      "    value.                                                         \n"
      "    Ex: \"max\" : 0.1  ( -->  [0.1, 0.1, ..., 0.1])                \n"
      "    If array, dimension must be equal to the \"axes\" dimension and\n"
      "    `\"sym_axes\" must be false`.                                  \n"
      "    Ex: \"max\" : [0.05, 0.1, 0.1]                                 \n\n"

      "  increment: number, or array of numbers (optional)                \n"
      "    Amount by which to increment counter elements. If number,      \n"
      "    specifies using a constant array of DoF space dimension with   \n"
      "    that given value.                                              \n"
      "    Ex: \"increment\" : 0.01  ( -->  [0.01, 0.01, ..., 0.01])      \n"
      "    If array, dimension must be equal to the \"axes\" dimension and\n"
      "    `\"sym_axes\" must be false`.                                  \n"
      "    Ex: \"increment\" : [0.005, 0.01, 0.01]                        \n"
      "    One of \"increment\" or \"num\" must be given.                 \n\n"

      "  num: int, or array of int (optional)                             \n"
      "    Number of values to include. Must be >= 1. If \"num\" is 1,    \n"
      "    only include a point at the \"min\" value along specified      \n"
      "    dimensions (this is equivalent to min=min, increment=(max-min),\n"
      "    max=max-increment/10.). If \"num\" is >1, include that many    \n"
      "    points, including the \"min\" and \"max\" (this is equivalent  \n"
      "    to min=min, increment=(max-min)/(num-1), max=max+inc/10.).     \n"
      "    Ex: \"increment\" : 11  ( --> [11, 11, ..., 11])               \n"
      "    If array, dimension must be equal to the \"axes\" dimension.   \n"
      "    Ex: \"num\" : [5, 11, 11]                                      \n"
      "    One of \"increment\" or \"num\" must be given.                 \n\n"

      "  min_nonzero: integer (optional, default = 0)                     \n"
      "    Minimum number of coordinate amplitudes that are allowed       \n"
      "    to be nonzero. Must be less than or equal to \"max_nonzero\".  \n\n"

      "  max_nonzero: integer (optional, default = axes.rows())           \n"
      "    Maximum number of coordinate amplitudes that are allowed to be \n"
      "    nonzero. Must be less than or equal to the \"axes\" dimension. \n\n"

      "  output_dir: string (optional, default=current path)              \n"
      "    Selects where output files are written.                        \n\n";

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

  return name() + ": \n\n" + description + custom_options +
         standard_ConfigEnumInput_help() + examples;
}

std::string ConfigEnumSiteDoFsInterface::name() const {
  return ConfigEnumSiteDoFs::enumerator_name;
}

/// Parse the ConfigEnumSiteDoFsParams JSON input
void parse(InputParser<ConfigEnumSiteDoFsParams> &parser,
           ConfigEnumInput const &initial_state, AxesCounterParams &axes_params,
           bool sym_axes_option, bool exclude_homogeneous_modes) {
  parser.value = notstd::make_unique<ConfigEnumSiteDoFsParams>();
  auto &params = *parser.value;

  // 1) get DoF type -------------------------------------------------
  // "dof" -> params.dof
  parser.require(params.dof, "dof");

  Supercell const &supercell = initial_state.configuration().supercell();
  Index dof_space_dimension = get_dof_space_dimension(
      params.dof, supercell.prim(),
      supercell.sym_info().transformation_matrix_to_super(),
      initial_state.sites());

  // 2) get axes and normal coordinate grid
  // ------------------------------------------
  auto grid_parser = parser.parse_as<AxesCounterParams>(dof_space_dimension);
  if (grid_parser->valid()) {
    axes_params = *grid_parser->value;
  }

  if ((sym_axes_option || exclude_homogeneous_modes) &&
      !axes_params.scalar_input) {
    std::stringstream msg;
    msg << "Error: Vector input for enumeration ranges (\"min\", \"max\", "
           "\"increment\" or \"num\") is not allowed with `\"sym_axes\": true` "
           "or `\"exclude_homogeneous_modes\": true`.";
    throw std::runtime_error(msg.str());
  }

  // 3) get min/max nonzero amplitudes -----------------------------------

  // "min_nonzero" -> params.min_nonzero
  parser.optional_else(params.min_nonzero, "min_nonzero", Index{0});

  // "max_nonzero" -> params.max_nonzero
  // note that help indicates default==axes.rows(), but that is
  // params.axes.cols()
  parser.optional_else(params.max_nonzero, "max_nonzero",
                       Index{params.axes.cols()});

  // Throw error if min_nonzero exceeds max_nonzero
  if (params.min_nonzero > params.max_nonzero) {
    std::stringstream msg;
    msg << "'min_nonzero' value exceeds 'max_nonzero' value" << std::endl;
    parser.error.insert(msg.str());
  }
}

void require_all_input_have_the_same_number_of_selected_sites(
    InputParser<std::vector<std::pair<std::string, ConfigEnumInput>>> &parser) {
  // TODO: should this check if the supercells & sites are the same, and not
  // just the same number? if sym_axes==false, then using different sites /
  // supercells might be bad

  if (parser.valid() && parser.value != nullptr) {
    auto const &named_initial_states = *parser.value;
    Index nsites = named_initial_states[0].second.sites().size();
    for (auto const &name_value_pair : named_initial_states) {
      if (name_value_pair.second.sites().size() != nsites) {
        std::stringstream msg;
        msg << "Error in ConfigEnumSiteDoFs: Starting configurations or "
               "supercells passed to "
               "must all have the same number of selected sites.";
        parser.error.insert(msg.str());
      }
    }
  }
}

namespace ConfigEnumSiteDoFsInterface_impl {

typedef ConfigEnumData<ConfigEnumSiteDoFs, ConfigEnumInput> ConfigEnumDataType;

// This functor constructs a ConfigEnumSiteDoFs enumerator for a given initial
// state
struct MakeEnumerator {
  MakeEnumerator(ConfigEnumOptions const &_options,
                 ConfigEnumSiteDoFsParams const &_params,
                 AxesCounterParams const &_axes_params,
                 bool _make_symmetry_adapted_axes,
                 bool _exclude_homogeneous_modes,
                 DoFSpaceIO::SequentialDirectoryOutput &_dof_space_output)
      : log(CASM::log()),
        options(_options),
        params_template(_params),
        axes_params(_axes_params),
        make_symmetry_adapted_axes(_make_symmetry_adapted_axes),
        exclude_homogeneous_modes(_exclude_homogeneous_modes),
        calc_wedges(false),
        dof_space_output(_dof_space_output) {}

  Log &log;
  ConfigEnumOptions const &options;
  ConfigEnumSiteDoFsParams const &params_template;
  AxesCounterParams const &axes_params;
  bool make_symmetry_adapted_axes;
  bool exclude_homogeneous_modes;
  bool calc_wedges;
  DoFSpaceIO::SequentialDirectoryOutput &dof_space_output;

  // constructs a ConfigEnumSiteDoFs for each initial_state
  ConfigEnumSiteDoFs operator()(Index index, std::string name,
                                ConfigEnumInput const &initial_state) const;

  // constructs a symmetry adapted dof space and writes it to file as a record
  DoFSpace make_and_write_dof_space(
      Index index, std::string name, ConfigEnumInput const &initial_state,
      std::optional<VectorSpaceSymReport> &sym_report) const;

  // constructs a DataFormatter to record enumeration results
  DataFormatter<ConfigEnumDataType> make_formatter() const;
};

// constructs a ConfigEnumSiteDoFs for each initial_state
ConfigEnumSiteDoFs MakeEnumerator::operator()(
    Index index, std::string name, ConfigEnumInput const &initial_state) const {
  std::optional<VectorSpaceSymReport> sym_report;
  DoFSpace dof_space =
      make_and_write_dof_space(index, name, initial_state, sym_report);

  ConfigEnumSiteDoFsParams params = params_template;
  params.axes = dof_space.basis();

  // set enumeration ranges
  if (axes_params.scalar_input) {
    int dim = dof_space.dim();
    params.min_val = Eigen::VectorXd::Constant(dim, axes_params.min_scalar);
    params.max_val = Eigen::VectorXd::Constant(dim, axes_params.max_scalar);
    params.inc_val = Eigen::VectorXd::Constant(dim, axes_params.inc_scalar);
  } else {
    params.min_val = axes_params.min_vector;
    params.max_val = axes_params.max_vector;
    params.inc_val = axes_params.inc_vector;
  }

  return ConfigEnumSiteDoFs{initial_state, params};
}

// constructs a symmetry adapted dof space and writes it to file as a record
DoFSpace MakeEnumerator::make_and_write_dof_space(
    Index index, std::string name, ConfigEnumInput const &initial_state,
    std::optional<VectorSpaceSymReport> &sym_report) const {
  DoFSpace dof_space =
      make_dof_space(params_template.dof, initial_state, axes_params.axes);
  if (make_symmetry_adapted_axes) {
    log << "Performing DoF space analysis: " << name << std::endl;
    auto const &sym_info = initial_state.configuration().supercell().sym_info();
    std::vector<PermuteIterator> group = make_invariant_subgroup(initial_state);
    if (exclude_homogeneous_modes) {
      dof_space = exclude_homogeneous_mode_space(dof_space);
    }
    dof_space_output.write_symmetry(index, name, initial_state, group);
    dof_space = make_symmetry_adapted_dof_space(dof_space, sym_info, group,
                                                calc_wedges, sym_report);
  }
  dof_space_output.write_dof_space(index, dof_space, name, initial_state,
                                   sym_report);
  return dof_space;
}

// constructs a DataFormatter to record enumeration results
DataFormatter<ConfigEnumDataType> MakeEnumerator::make_formatter() const {
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

  return formatter;
}

}  // namespace ConfigEnumSiteDoFsInterface_impl

void ConfigEnumSiteDoFsInterface::run(
    PrimClex &primclex, jsonParser const &json_options,
    jsonParser const &cli_options_as_json) const {
  Log &log = CASM::log();

  log.subsection().begin("ConfigEnumSiteDoFs");
  ParentInputParser parser =
      make_enum_parent_parser(log, json_options, cli_options_as_json);
  std::runtime_error error_if_invalid{
      "Error reading ConfigEnumSiteDoFs JSON input"};

  log.custom("Checking input");

  // 1) Parse ConfigEnumOptions ------------------

  auto options_parser_ptr = parser.parse_as<ConfigEnumOptions>(
      ConfigEnumSiteDoFs::enumerator_name, primclex,
      primclex.settings().query_handler<Configuration>().dict());
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  ConfigEnumOptions const &options = *options_parser_ptr->value;
  print_options(log, options);
  log.set_verbosity(options.verbosity);

  // 2) Parse initial enumeration states ------------------
  typedef std::vector<std::pair<std::string, ConfigEnumInput>>
      NamedInitialEnumerationStates;
  auto input_parser_ptr = parser.parse_as<NamedInitialEnumerationStates>(
      primclex.shared_prim(), &primclex, primclex.db<Supercell>(),
      primclex.db<Configuration>());
  require_all_input_have_the_same_number_of_selected_sites(*input_parser_ptr);
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  auto const &named_initial_states = *input_parser_ptr->value;
  print_initial_states(log, named_initial_states);

  // 3) Parse ConfigEnumSiteDoFsParams ------------------

  // 3a) check for "sym_axes" option:
  bool sym_axes_option;
  parser.optional_else(sym_axes_option, "sym_axes", false);
  log.indent() << "sym_axes: " << std::boolalpha << sym_axes_option
               << std::endl;

  // 3c) check for "exclude_homogeneous_modes" option
  bool exclude_homogeneous_modes;
  parser.optional_else(exclude_homogeneous_modes, "exclude_homogeneous_modes",
                       false);
  log.indent() << "exclude_homogeneous_modes: " << std::boolalpha
               << exclude_homogeneous_modes << std::endl;

  // axes: column vector matrix
  // - user input, else default = identity matrix of initial state dof space
  //   dimension
  // - If sym_axes==false: normal mode coordinates used directly
  // - If sym_axes==true: subspace to be symmetrized
  AxesCounterParams axes_params;
  auto params_parser_ptr = parser.parse_as<ConfigEnumSiteDoFsParams>(
      named_initial_states[0].second, axes_params, sym_axes_option,
      exclude_homogeneous_modes);
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  ConfigEnumSiteDoFsParams const &params = *params_parser_ptr->value;
  log.indent() << "axes: (column vectors) \n" << axes_params.axes << std::endl;
  if (axes_params.scalar_input) {
    log.indent() << "min: " << axes_params.min_scalar << std::endl;
    log.indent() << "max: " << axes_params.max_scalar << std::endl;
    log.indent() << "increment: " << axes_params.inc_scalar << std::endl;
  } else {
    log.indent() << "min: " << axes_params.min_vector.transpose() << std::endl;
    log.indent() << "max: " << axes_params.max_vector.transpose() << std::endl;
    log.indent() << "increment: " << axes_params.inc_vector.transpose()
                 << std::endl;
  }
  log.indent() << "min_nonzero: " << params.min_nonzero << std::endl;
  log.indent() << "max_nonzero: " << params.max_nonzero << std::endl;

  // 3b) check for "print_dof_space_and_quit" option:
  bool print_dof_space_and_quit_option;
  parser.optional_else(print_dof_space_and_quit_option,
                       "print_dof_space_and_quit", false);

  // parse "output_dir" (optional, default = current_path)
  fs::path output_dir;
  parser.optional_else(output_dir, "output_dir", fs::current_path());
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  log << std::endl;

  // 4) Enumerate configurations ------------------

  DoFSpaceIO::SequentialDirectoryOutput dof_space_output{output_dir};
  ConfigEnumSiteDoFsInterface_impl::MakeEnumerator make_enumerator_f{
      options,
      params,
      axes_params,
      sym_axes_option,
      exclude_homogeneous_modes,
      dof_space_output};

  if (print_dof_space_and_quit_option) {
    log.begin("Print DoF Space and Quit Option");
    log << "For large spaces this may be slow..." << std::endl;
    std::optional<VectorSpaceSymReport> sym_report;
    Index i = 0;
    for (auto const &pair : named_initial_states) {
      std::string const &name = pair.first;
      ConfigEnumInput const &initial_state = pair.second;
      make_enumerator_f.make_and_write_dof_space(i, name, initial_state,
                                                 sym_report);
      i++;
    }
    log.end_section();
    return;
  }

  log << std::endl;
  log.begin("ConfigEnumSiteDoFs enumeration");
  enumerate_configurations(
      primclex, options, make_enumerator_f, named_initial_states.begin(),
      named_initial_states.end(), make_enumerator_f.make_formatter());
  log.end_section();
}

}  // namespace CASM
