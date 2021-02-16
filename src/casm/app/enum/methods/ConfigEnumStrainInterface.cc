#include "casm/app/enum/methods/ConfigEnumStrainInterface.hh"

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
#include "casm/clex/ConfigEnumStrain.hh"
#include "casm/clex/ConfigIOStrain.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/DoFSpace.hh"
#include "casm/enumerator/io/dof_space_analysis.hh"
#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"
#include "casm/symmetry/SymRepTools.hh"
#include "casm/symmetry/io/json/SymRepTools.hh"

namespace CASM {

// namespace adapter {
//
//   Configuration const &
//   Adapter<Configuration, ConfigEnumData<ConfigEnumStrain,
//   ConfigEnumInput>>::operator()(
//     ConfigEnumData<ConfigEnumStrain, ConfigEnumInput> const &adaptable) const
//     { return adaptable.configuration;
//   }
//
//   Supercell const &
//   Adapter<Supercell, ConfigEnumData<ConfigEnumStrain,
//   ConfigEnumInput>>::operator()(
//     ConfigEnumData<ConfigEnumStrain, ConfigEnumInput> const &adaptable) const
//     { return adaptable.configuration.supercell();
//   }
//
// }

namespace ConfigEnumIO {
/// Template specialization to get current normal coordinate from enumerator.
template <>
Eigen::VectorXd get_normal_coordinate(ConfigEnumStrain const &enumerator) {
  return enumerator.normal_coordinate();
}

template <typename ConfigEnumDataType>
GenericDatumFormatter<Index, ConfigEnumDataType> subwedge_index() {
  return GenericDatumFormatter<Index, ConfigEnumDataType>(
      "subwedge_index",
      "Subwedge index (determines axes defining normal coordinate).",
      [](ConfigEnumDataType const &data) {
        return data.enumerator.subwedge_index();
      });
}
}  // namespace ConfigEnumIO

std::string ConfigEnumStrainInterface::desc() const {
  std::string description =
      "Generate strain perturbations of one or more initial "
      "configurations.\n\n";

  std::string custom_options =
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

      "  print_dof_space_and_quit: boolean (optional, default=false)      \n"
      "    If true, print DoF space information for each initial          \n"
      "    enumeration state and quit. If `\"sym_axes\": true`, will also \n"
      "    print irreducible subspaces and symmetry-adapted axes.         \n\n"

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

      "  trim_corners: bool (optional, default=true)                      \n"
      "    If true, any grid points outside the largest ellipsoid         \n"
      "    inscribed within the extrema of the grid will be discarded.    \n\n"

      "  output_dir: string (optional, default=current path)              \n"
      "    Selects where output files are written.                        \n\n";

  std::string examples =
      "  Examples:\n"
      "    To enumerate all strain perturbations of a particular "
      "configuration:\n"
      "      casm enum --method ConfigEnumStrain -i \n"
      "      '{ \n"
      "         \"confignames\" : [\"SCEL4_1_4_1_0_0_0/3\"],\n"
      "         \"increment\" : [0.01, 0.01, 0.01, 0., 0., 0.],\n"
      "         \"max\" : [0.05, 0.05, 0,05, 0., 0., 0.]\n"
      "       }' \n\n";

  return name() + ": \n\n" + description + custom_options +
         standard_ConfigEnumInput_help() + examples;
}

std::string ConfigEnumStrainInterface::name() const {
  return ConfigEnumStrain::enumerator_name;
}

namespace ConfigEnumStrainInterface_impl {

typedef ConfigEnumData<ConfigEnumStrain, ConfigEnumInput> ConfigEnumDataType;

// Holds parsed input parameters and acts as functor to construct
// a ConfigEnumStrain enumerator for each initial state
struct MakeEnumerator {
  MakeEnumerator(ConfigEnumOptions const &_options,
                 ConfigEnumStrainParams const &_params,
                 AxesCounterParams const &_axes_params,
                 bool _make_symmetry_adapted_axes,
                 DoFSpaceIO::SequentialDirectoryOutput &_dof_space_output)
      : log(CASM::log()),
        options(_options),
        params_template(_params),
        axes_params(_axes_params),
        make_symmetry_adapted_axes(_make_symmetry_adapted_axes),
        calc_wedges(true),
        dof_space_output(_dof_space_output) {}

  Log &log;
  ConfigEnumOptions const &options;
  ConfigEnumStrainParams const &params_template;
  AxesCounterParams const &axes_params;
  bool make_symmetry_adapted_axes;
  bool calc_wedges;
  DoFSpaceIO::SequentialDirectoryOutput &dof_space_output;

  // constructs a ConfigEnumStrain for each initial_state
  ConfigEnumStrain operator()(Index index, std::string name,
                              ConfigEnumInput const &initial_state) const;

  // constructs a symmetry adapted dof space and writes it to file as a record
  DoFSpace make_and_write_dof_space(
      Index index, std::string name, ConfigEnumInput const &initial_state,
      std::optional<VectorSpaceSymReport> &sym_report) const;

  // constructs a DataFormatter to record enumeration results
  DataFormatter<ConfigEnumDataType> make_formatter() const;
};

// constructs a ConfigEnumStrain for each initial_state
ConfigEnumStrain MakeEnumerator::operator()(
    Index index, std::string name, ConfigEnumInput const &initial_state) const {
  std::optional<VectorSpaceSymReport> sym_report;
  DoFSpace dof_space =
      make_and_write_dof_space(index, name, initial_state, sym_report);

  ConfigEnumStrainParams params = params_template;
  if (make_symmetry_adapted_axes) {  // if sym_axes==true, use wedges
    params.wedges = sym_report->irreducible_wedge;
  } else {
    params.wedges.clear();
    params.wedges.push_back(
        SymRepTools::SubWedge::make_dummy(dof_space.basis()));
  }

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

  return ConfigEnumStrain{initial_state, params};
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
  std::string prim_strain_metric = xtal::get_strain_metric(params_template.dof);
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
      ConfigEnumIO::initial_state_configname<ConfigEnumDataType>());
  if (make_symmetry_adapted_axes) {
    formatter.push_back(ConfigEnumIO::subwedge_index<ConfigEnumDataType>());
  }
  formatter.push_back(
      ConfigEnumIO::normal_coordinate<ConfigEnumDataType>(),
      make_datum_formatter_adapter<ConfigEnumDataType, Configuration>(
          ConfigIO::DoFStrain(prim_strain_metric)));
  if (prim_strain_metric != "F") {  // may not be necessary
    formatter.push_back(
        make_datum_formatter_adapter<ConfigEnumDataType, Configuration>(
            ConfigIO::DoFStrain("F")));
  }
  return formatter;
}

}  // namespace ConfigEnumStrainInterface_impl

/// Parse the ConfigEnumStrainParams (except parse "axes" instead of wedges)
/// from JSON input
///
/// \param parser InputParser will contain ConfigEnumStrainParams if successful
/// and error messages
///               otherwise
/// \param axes Eigen::MatriXd speciying a custom user subspace / normal mode
/// axes if provied,
///        otherwise set to Identity matrix of size dof_space_dimension x
///        dof_space_dimension.
/// \param dof_space_dimension DoF space dimension detected from initial
/// enumeration states
void parse(InputParser<ConfigEnumStrainParams> &parser,
           ConfigEnumInput const &initial_state, AxesCounterParams &axes_params,
           bool sym_axes_option) {
  parser.value = notstd::make_unique<ConfigEnumStrainParams>();
  auto &params = *parser.value;

  Supercell const &supercell = initial_state.configuration().supercell();

  // 1) get DoF type -------------------------------------------------
  // "dof" -> params.dof
  try {
    params.dof = xtal::get_strain_dof_key(supercell.prim());
  } catch (std::exception &e) {
    parser.error.insert(e.what());
    return;
  }

  Index dof_space_dimension = get_dof_space_dimension(
      params.dof, supercell.prim(),
      supercell.sym_info().transformation_matrix_to_super(),
      initial_state.sites());

  // 2) get axes and enumeration counter parameters -------------------
  auto grid_parser = parser.parse_as<AxesCounterParams>(dof_space_dimension);
  if (grid_parser->valid()) {
    axes_params = *grid_parser->value;
  }

  if (sym_axes_option && !axes_params.scalar_input) {
    std::stringstream msg;
    msg << "Error: Vector input for enumeration ranges (\"min\", \"max\", "
           "\"increment\" or \"num\") is not allowed with `\"sym_axes\": "
           "true`.";
    throw std::runtime_error(msg.str());
  }

  // 3) set auto_range option  ---------------------------------------
  if (!parser.self.contains("min") && sym_axes_option == true) {
    params.auto_range = true;
  } else {
    params.auto_range = false;
  }

  // 4) get trim_corners option  -------------------------------------
  // "trim_corners" -> params.trim_corners (bool, default true)
  parser.optional_else(params.trim_corners, "trim_corners", true);
}

void ConfigEnumStrainInterface::run(
    PrimClex &primclex, jsonParser const &json_options,
    jsonParser const &cli_options_as_json) const {
  using namespace ConfigEnumStrainInterface_impl;

  Log &log = CASM::log();

  log.subsection().begin("ConfigEnumStrain");
  ParentInputParser parser =
      make_enum_parent_parser(log, json_options, cli_options_as_json);
  std::runtime_error error_if_invalid{
      "Error reading ConfigEnumStrain JSON input"};

  log.custom("Checking input");

  // 1) Parse ConfigEnumOptions ------------------
  auto options_parser_ptr = parser.parse_as<ConfigEnumOptions>(
      ConfigEnumStrain::enumerator_name, primclex,
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
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  auto const &named_initial_states = *input_parser_ptr->value;
  print_initial_states(log, named_initial_states);

  // 3) Parse ConfigEnumStrainParams ------------------

  // 3a) check for "sym_axes" option:
  bool sym_axes_option;
  parser.optional_else(sym_axes_option, "sym_axes", false);
  log.indent() << "sym_axes: " << std::boolalpha << sym_axes_option
               << std::endl;

  // 3b) parse ConfigEnumStrainParams (except parse "axes" instead of "wedges")

  // axes: column vector matrix
  // - user input, else default = identity matrix of initial state dof space
  //   dimension
  // - If sym_axes==false: normal mode coordinates used directly
  // - If sym_axes==true: subspace to be symmetrized
  AxesCounterParams axes_params;
  auto params_parser_ptr = parser.parse_as<ConfigEnumStrainParams>(
      named_initial_states[0].second, axes_params, sym_axes_option);
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  ConfigEnumStrainParams const &params = *params_parser_ptr->value;
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
  log.indent() << "auto_range: " << std::boolalpha << params.auto_range
               << std::endl;
  log.indent() << "trim_corners: " << std::boolalpha << params.trim_corners
               << std::endl;

  // 3c) check for "print_dof_space_and_quit" option:
  bool print_dof_space_and_quit_option;
  parser.optional_else(print_dof_space_and_quit_option,
                       "print_dof_space_and_quit", false);
  log.indent() << "print_dof_space_and_quit: " << std::boolalpha
               << print_dof_space_and_quit_option << std::endl;

  // 3d) parse "output_dir" (optional, default = current_path)
  fs::path output_dir;
  parser.optional_else(output_dir, "output_dir", fs::current_path());
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  log << std::endl;

  // 4) Enumerate configurations ------------------
  DoFSpaceIO::SequentialDirectoryOutput dof_space_output{output_dir};
  MakeEnumerator make_enumerator_f{options, params, axes_params,
                                   sym_axes_option, dof_space_output};

  if (print_dof_space_and_quit_option) {
    log.begin("Print DoF Space and Quit Option");
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
  log.begin("ConfigEnumStrain enumeration");
  enumerate_configurations(
      primclex, options, make_enumerator_f, named_initial_states.begin(),
      named_initial_states.end(), make_enumerator_f.make_formatter());
  log.end_section();
}

}  // namespace CASM
