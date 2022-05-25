#include "casm/enumerator/DoFSpace.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/optional.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/io/json/ConfigDoF_json_io.hh"
#include "casm/clex/io/json/Configuration_json_io.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/io/SimpleStructureIO.hh"
#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/symmetry/SupercellSymInfo.hh"
#include "casm/symmetry/io/json/SymRepTools.hh"

namespace CASM {

namespace {  // implementation for to_json(DoFSpace const &dofspace, ...)

/// Note: dofspace must be local continuous DoFSpace, with
///     `dofspace.include_all_sites()==true, else will throw
void _add_homogeneous_mode_info(jsonParser &json, DoFSpace const &dofspace) {
  jsonParser &irreps_json = json["irreducible_representations"];

  Eigen::MatrixXd homogeneous_mode_space =
      make_homogeneous_mode_space(dofspace);

  auto print_string = [&](std::vector<Index> const &indices) {
    std::vector<std::string> string_of_indices;
    for (auto &i : indices) {
      std::string s = "q" + to_sequential_string(i + 1, dofspace.dim());
      string_of_indices.push_back(s);
    }
    return string_of_indices;
  };

  VectorSpaceMixingInfo mixing_info(dofspace.basis(), homogeneous_mode_space,
                                    CASM::TOL);
  irreps_json["adapted_axes_which_are_not_homogeneous_modes"] =
      print_string(mixing_info.axes_not_in_subspace);
  irreps_json["adapted_axes_mixed_with_homogeneous_modes"] =
      print_string(mixing_info.axes_mixed_with_subspace);
  irreps_json["adapted_axes_which_are_homogeneous_modes"] =
      print_string(mixing_info.axes_in_subspace);
}
}  // namespace

void from_json(DoFSpace &dofspace, jsonParser const &json,
               std::shared_ptr<Structure const> const &shared_prim) {
  dofspace = jsonConstructor<DoFSpace>::from_json(json, shared_prim);
}

jsonParser &to_json(DoFSpace const &dofspace, jsonParser &json,
                    std::optional<std::string> const &identifier,
                    std::optional<ConfigEnumInput> const &input_state,
                    std::optional<VectorSpaceSymReport> const &sym_report) {
  json["dof"] = dofspace.dof_key();
  json["transformation_matrix_to_supercell"] =
      dofspace.transformation_matrix_to_super();
  json["sites"] = dofspace.sites();
  json["basis"] = dofspace.basis().transpose();
  json["glossary"] = dofspace.axis_glossary();
  json["axis_site_index"] = dofspace.axis_site_index();
  json["axis_dof_component"] = dofspace.axis_dof_component();

  if (identifier.has_value()) {
    json["identifier"] = identifier;
  }
  if (input_state.has_value()) {
    json["state"] = input_state;
  }
  if (sym_report.has_value()) {
    to_json(sym_report, json);
    if (dofspace.includes_all_sites() && dofspace.dof_key() == "disp") {
      _add_homogeneous_mode_info(json, dofspace);
    }
  }
  return json;
}

jsonParser &to_json(
    DoFSpace const &dofspace, jsonParser &json,
    std::optional<std::string> const &identifier,
    std::optional<ConfigEnumInput> const &input_state,
    std::optional<SymRepTools_v2::VectorSpaceSymReport> const &sym_report) {
  json["dof"] = dofspace.dof_key();
  json["transformation_matrix_to_supercell"] =
      dofspace.transformation_matrix_to_super();
  json["sites"] = dofspace.sites();
  json["basis"] = dofspace.basis().transpose();
  json["glossary"] = dofspace.axis_glossary();
  json["axis_site_index"] = dofspace.axis_site_index();
  json["axis_dof_component"] = dofspace.axis_dof_component();

  if (identifier.has_value()) {
    json["identifier"] = identifier;
  }
  if (input_state.has_value()) {
    json["state"] = input_state;
  }
  if (sym_report.has_value()) {
    to_json(sym_report, json);
    if (dofspace.includes_all_sites() && dofspace.dof_key() == "disp") {
      _add_homogeneous_mode_info(json, dofspace);
    }
  }
  return json;
}

DoFSpace jsonConstructor<DoFSpace>::from_json(
    jsonParser const &json,
    std::shared_ptr<Structure const> const &shared_prim) {
  InputParser<DoFSpace> parser{json, shared_prim};

  std::runtime_error error_if_invalid{"Error reading DoFSpace from JSON"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

  return std::move(*parser.value);
}

std::unique_ptr<DoFSpace> jsonMake<DoFSpace>::make_from_json(
    jsonParser const &json,
    std::shared_ptr<Structure const> const &shared_prim) {
  InputParser<DoFSpace> parser{json, shared_prim};

  std::runtime_error error_if_invalid{"Error reading DoFSpace from JSON"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

  return std::move(parser.value);
}

void parse(InputParser<DoFSpace> &parser,
           std::shared_ptr<Structure const> const &shared_prim) {
  DoFKey dof_key;
  parser.require(dof_key, "dof");

  std::optional<Eigen::Matrix3l> transformation_matrix_to_super;
  parser.optional(transformation_matrix_to_super,
                  "transformation_matrix_to_supercell");

  std::optional<std::set<Index>> sites;
  parser.optional(sites, "sites");

  std::optional<Eigen::MatrixXd> basis;
  parser.optional(basis, "basis");
  if (basis.has_value()) {
    basis = basis.value().transpose();
  }

  if (parser.valid()) {
    parser.value = notstd::make_unique<DoFSpace>(
        shared_prim, dof_key, transformation_matrix_to_super, sites, basis);
  }
}

namespace {  // implementation for parse_axes_counter_range

/// For 'num input, updates 'inc' and 'max', based on 'min, 'inc', and 'num'
///
/// If 'num' < 1: insert error
/// If 'num' == 1: min=min; inc=max-min; max=max-inc/10.;
/// If 'num' > 1: min=min; inc=(max-min)/(num-1); max=max+inc/10.;
void _enforce_num(InputParser<AxesCounterParams> &parser, double &inc,
                  double &max, double min, int num) {
  if (num < 1) {
    parser.insert_error("num", "Error: 'num' must be >= 1.");
  } else if (num == 1) {
    inc = (max - min);
    max -= inc / 10.;
  } else {
    inc = (max - min) / (num - 1);
    max += inc / 10.;
  }
}

/// Parse scalar "min", "max", "increment"/"num"
void _parse_scalar_input(InputParser<AxesCounterParams> &parser) {
  AxesCounterParams &params = *parser.value;
  params.scalar_input = true;

  parser.optional_else(params.min_scalar, "min", 0.);
  parser.require(params.max_scalar, "max");

  if (parser.self.contains("increment") == parser.self.contains("num")) {
    parser.error.insert("Error: require one of 'increment' or 'num'.");
    return;
  }
  if (parser.self.contains("increment")) {
    parser.optional(params.inc_scalar, "increment");
  } else {
    /// Alternative to "increment", specify number of values along dimension
    /// (including min and max)
    int num_scalar;
    parser.optional(num_scalar, "num");
    _enforce_num(parser, params.inc_scalar, params.max_scalar,
                 params.min_scalar, num_scalar);
  }
}

/// Parse vector "min", "max", "increment"/"num"
void _parse_vector_input(InputParser<AxesCounterParams> &parser) {
  AxesCounterParams &params = *parser.value;
  int dim = params.axes.cols();
  params.scalar_input = false;

  // "min"
  Eigen::VectorXd default_min = Eigen::VectorXd::Constant(dim, 0.);
  parser.optional_else(params.min_vector, "min", default_min);
  if (params.min_vector.size() != dim) {
    parser.insert_error("min", "Error: 'min' size != axes dimension");
  }

  // "max"
  parser.require(params.max_vector, "max");
  if (params.max_vector.size() != dim) {
    parser.insert_error("max", "Error: 'max' size != axes dimension");
  }

  // "increment" or "num"
  if (parser.self.contains("increment") == parser.self.contains("num")) {
    parser.error.insert("Error: require one of 'increment' or 'num'.");
    return;
  }
  if (parser.self.contains("increment")) {
    parser.optional(params.inc_vector, "increment");
    if (params.inc_vector.size() != dim) {
      parser.insert_error("increment",
                          "Error: 'increment' size != axes dimension");
    }
  } else {
    /// Alternative to "increment", specify number of values along dimension
    /// (including min and max)
    Eigen::VectorXi num_vector;
    parser.optional(num_vector, "num");
    if (num_vector.size() != dim) {
      parser.insert_error("num", "Error: 'num' size != axes dimension");
    }
    for (int i = 0; i < num_vector.size(); ++i) {
      _enforce_num(parser, params.inc_vector[i], params.max_vector[i],
                   params.min_vector[i], num_vector[i]);
    }
  }
}

}  // namespace

/// Parse "min", "max", and one of "increment" or "num"
///
/// \param parser AxesCounterParams parser
///
/// Examples of equivalent input for a grid in 2d space:
/// - `{ "min": 0., "max": 0.101, "increment": 0.01}`
/// - `{ "min": 0., "max": 0.1, "num": 11}`
/// - `{ "min": [0., 0.], "max": [0.101, 0.101], "increment": [0.1, 0.1]}`
/// - `{ "min": [0., 0.], "max": [0.1, 0.1], "num": [11, 11]}`
///
/// Notes:
/// - all range inputs must be consistently either scalar or array
/// - "min", "max", "increment" must be number or array of number
/// - "num" must be int or array of int, specifies number of values to count
///   over (including min and max values, unless num==1, then exclude max)
/// - only "min" is optional, with default value scalar 0, or vector 0,
///   depending on type of "max"
/// - only one of "increment", "num" may be present
/// - sets AxesCounterParams members `scalar_input` based on input
void parse_axes_counter_range(InputParser<AxesCounterParams> &parser) {
  if (parser.value == nullptr) {
    throw std::runtime_error("Unknown AxesCounterParams parsing error");
    return;
  }

  // Use "max" to check for scalar or array input
  if (!parser.self.contains("max")) {
    parser.insert_error("max", "Error: missing required parameter 'max'.");
    return;
  } else if (parser.self["max"].is_number()) {
    _parse_scalar_input(parser);
  } else if (parser.self["max"].is_array()) {
    _parse_vector_input(parser);
  } else {
    std::stringstream msg;
    msg << "Error: Parameter 'max' must be a number, or an "
           "array of numbers matching the 'axes' dimension.";
    parser.insert_error("max", msg.str());
    return;
  }
}

/// Parse DoF space axes from a JSON object
///
/// Expected format: (this parses "axes" only)
///
/// Full space case:
/// {
///   "axes": {
///     "q1": [q10, q11, q12, ...], /// note: starts with "q1", not "q0"
///     "q2": [q20, q21, q22, ...],
///     "q3": [q30, q31, q32, ...],
///     ...
///   },
///   "min": [q1min, q2min, q3min, ...], /// optional, default is zeros vector
///   "max": [q1max, q2max, q3max, ...], /// required
///   "increment": [q1inc, q2inc, q3inc, ...],  /// required
/// }
///
/// Subspace case:
/// - if some qi (i in range [1, dof_space_dimension]) are missing, then use the
///   subspace specified by the axes that are provided
/// - if scalars, min, max, increment will be used along each dimension
/// - if arrays, min.size() / max.size() / increment.size() must equal
///   axes.size()
/// - example:
/// {
///   "axes": {
///     "q1": [q10, q11, q12, ...],
///     "q2": [q20, q21, q22, ...],
///     "q3": [q30, q31, q32, ...],
///     "q5": [q50, q51, q52, ...],
///   },
///   "min": [q1min, q2min, q3min, q5min],        /// optional, default is zeros
///   vector "max": [q1max, q2max, q3max, q5max], /// required
///   "increment": [q1inc, q2inc, q3inc, q5inc],  /// required
/// }
///
/// notes:
/// - if total dof_space_dimension is double digit, then use "q01", "q02", ...
/// etc.
/// - if total dof_space_dimension is triple digit, then use "q001", "q002", ...
/// etc.
///
void parse_axes_from_object(InputParser<AxesCounterParams> &parser,
                            Index dof_space_dimension) {
  if (parser.value == nullptr) {
    throw std::runtime_error("Unknown AxesCounterParams parsing error");
    return;
  }
  Eigen::MatrixXd inaxes =
      Eigen::MatrixXd::Zero(dof_space_dimension, dof_space_dimension);

  std::set<Index> found;
  for (Index i = 0; i < dof_space_dimension; ++i) {
    std::string axis_name =
        "q" + to_sequential_string(i + 1, dof_space_dimension);
    auto subparser =
        parser.subparse_if<Eigen::VectorXd>(fs::path{"axes"} / axis_name);
    if (subparser->value != nullptr) {
      if (subparser->value->size() != dof_space_dimension) {
        std::stringstream msg;
        msg << "Error reading axis vector '" << axis_name
            << "': expected size=" << dof_space_dimension
            << " found size=" << subparser->value->size();
        subparser->error.insert(msg.str());
      } else {
        inaxes.col(found.size()) = *subparser->value;
      }
      found.insert(i);
    }
  }

  Index subspace_dimension = found.size();
  parser.value->axes = inaxes.leftCols(subspace_dimension);
}

/// Read "axes" from a row-vector JSON matrix, store in `parser.value->axes` as
/// a column vector matrix
void parse_axes_from_array(InputParser<AxesCounterParams> &parser,
                           Index dof_space_dimension) {
  if (parser.value == nullptr) {
    throw std::runtime_error("Unknown AxesCounterParams parsing error");
    return;
  }

  auto subparser = parser.subparse_if<Eigen::MatrixXd>("axes");
  if (!subparser->valid()) {
    return;
  }
  Eigen::MatrixXd &axes = parser.value->axes;
  axes = subparser->value->transpose();

  // check axes dimensions:
  if (axes.rows() != dof_space_dimension) {
    // Note: message "columns" refers to JSON input axes, transpose of
    // params.axes
    std::stringstream msg;
    msg << "Number of columns of 'axes' must be equal to site DoF space "
           "dimension ("
        << dof_space_dimension << "). Size as parsed: " << axes.rows();
    subparser->error.insert(msg.str());
  }
  if (axes.cols() > dof_space_dimension) {
    // Note: message "rows" refers to JSON input axes, transpose of params.axes
    std::stringstream msg;
    msg << "Number of coordinate axes (number of rows of 'axes') must be less "
           "than or equal to "
           "site DoF space dimension ("
        << dof_space_dimension << "). Number of axes parsed: " << axes.cols();
    subparser->error.insert(msg.str());
  }
}

/// Parse DoF space "axes"
///
/// Accepts a JSON object of axes or a row vector matrix of axes. May be rank
/// deficient.
///
/// \param parser JSON to be parsed
/// \param axes Set to be column vector matrix of axes
/// \param dof_space_dimension Full DoF space dimension. If no `"axes"` in the
///        JSON, the default axes are identity matrix of size
///        dof_space_dimension. Errors are inserted into parser if axes vectors
///        are not of length equal to dof_space_dimension.
///
/// TODO: document "axes" here
///
void parse_dof_space_axes(InputParser<AxesCounterParams> &parser,
                          Index dof_space_dimension) {
  if (parser.value == nullptr) {
    return;
  }
  Eigen::MatrixXd &axes = parser.value->axes;
  // "axes" -> axes
  if (!parser.self.contains("axes")) {
    axes = Eigen::MatrixXd::Identity(dof_space_dimension, dof_space_dimension);
  } else if (parser.self.contains("axes") && parser.self["axes"].is_obj()) {
    parse_axes_from_object(parser, dof_space_dimension);
  } else if (parser.self.contains("axes") && parser.self["axes"].is_array()) {
    parse_axes_from_array(parser, dof_space_dimension);
  } else {
    std::stringstream msg;
    msg << "The 'axes' must be a JSON object of axis vectors (named 'q1', "
           "'q2', etc.), or a row-vector matrix";
    parser.error.insert(msg.str());
  }
}

/// Parse DoFSpace subspace from  "axes" and normal coordinate grid counter from
/// "min", "max", and "increment"
void parse(InputParser<AxesCounterParams> &parser, Index dof_space_dimension) {
  parser.value = notstd::make_unique<AxesCounterParams>();
  parse_dof_space_axes(parser, dof_space_dimension);
  parse_axes_counter_range(parser);
  if (!parser.valid()) {
    parser.value.reset();
  }
}

}  // namespace CASM
