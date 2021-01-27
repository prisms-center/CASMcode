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
///
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
  parser.require(min, "min");
  parser.require(max, "max");
  parser.optional_else(dirs, "dirs", default_dirs);
  parser.optional_else(generating_matrix, "unit_cell", default_matrix);
  parser.value = notstd::make_unique<xtal::ScelEnumProps>(min, max + 1, dirs,
                                                          generating_matrix);
}
}  // namespace CASM
