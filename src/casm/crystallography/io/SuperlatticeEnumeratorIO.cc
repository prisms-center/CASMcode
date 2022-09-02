#include "casm/crystallography/io/SuperlatticeEnumeratorIO.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"

namespace CASM {
jsonParser &to_json(const xtal::ScelEnumProps &props, jsonParser &json) {
  json.put_obj();
  json["min"] = props.begin_volume();
  json["max"] = props.end_volume() - 1;
  json["dirs"] = props.dirs();
  json["unit_cell"] = props.generating_matrix();
  json["diagonal_only"] = props.diagonal_only();
  json["fixed_shape"] = props.fixed_shape();
  return json;
}

/// Make a ScelEnumProps object from JSON input
///
/// Expects a JSON object with options:
///     min: int (required)
///         Minimum volume supercells to enumerate
///     max: int (required)
///         Maximum volume supercells to enumerate
///     dirs: string (optional, default="abc")
///         Which lattice vectors of unit cell to enumerate over
///     unit_cell: 3x3 matrix of int (optional, default=identity matrix)
///         The unit cell to tile into supercells.
///     diagonal_only: bool (optional, default=false)
///         Restrict to diagonal multiples of unit_cell
///     fixed_shape: bool (optional, default=false)
///         Restrict to diagonal multiple of unit cell, with constant
///         cofficient for directions indicated by "dirs"
xtal::ScelEnumProps jsonConstructor<xtal::ScelEnumProps>::from_json(
    const jsonParser &json) {
  jsonParser tjson{json};
  InputParser<xtal::ScelEnumProps> parser{tjson};

  std::runtime_error error_if_invalid{
      "Error reading xtal::ScelEnumProps from JSON"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

  return *parser.value;
}

/// Make a ScelEnumProps object from JSON input
void from_json(xtal::ScelEnumProps &props, const jsonParser &json) {
  props = jsonConstructor<xtal::ScelEnumProps>::from_json(json);
}

/// Make a ScelEnumProps object from JSON input
void parse(InputParser<xtal::ScelEnumProps> &parser) {
  int min, max;
  std::string dirs;
  std::string default_dirs{"abc"};
  Eigen::Matrix3i generating_matrix;
  Eigen::Matrix3i default_matrix{Eigen::Matrix3i::Identity()};
  bool diagonal_only;
  bool fixed_shape;
  parser.require(min, "min");
  parser.require(max, "max");
  parser.optional_else(dirs, "dirs", default_dirs);
  parser.optional_else(generating_matrix, "unit_cell", default_matrix);
  parser.optional_else(diagonal_only, "diagonal_only", false);
  parser.optional_else(fixed_shape, "fixed_shape", false);
  parser.value = notstd::make_unique<xtal::ScelEnumProps>(
      min, max + 1, dirs, generating_matrix, diagonal_only, fixed_shape);
}
}  // namespace CASM
