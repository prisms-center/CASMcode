#include "casm/clex/io/json/Configuration_json_io.hh"

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/ConfigDoFTools.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/io/json/ConfigDoF_json_io.hh"

namespace CASM {

/// Read Configuration from JSON
///
/// See `from_json(Configuration &configuration, jsonParser const &json)` for
/// details.
Configuration jsonConstructor<Configuration>::from_json(
    jsonParser const &json,
    std::shared_ptr<const Structure> const &shared_prim) {
  return std::move(*jsonMake<Configuration>::make_from_json(json, shared_prim));
}

/// Read Configuration from JSON
///
/// See `from_json(Configuration &configuration, jsonParser const &json)` for
/// details.
std::unique_ptr<Configuration> jsonMake<Configuration>::make_from_json(
    jsonParser const &json,
    std::shared_ptr<const Structure> const &shared_prim) {
  auto &log = CASM::log();
  ParentInputParser parser{json};
  std::runtime_error error_if_invalid{
      "Error reading Configuration from JSON input"};

  Eigen::Matrix3l T;
  parser.require(T, "transformation_matrix_to_supercell");
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  auto shared_supercell = std::make_shared<Supercell const>(shared_prim, T);

  ConfigDoF configdof = make_configdof(*shared_supercell);
  parser.optional(configdof, "dof");

  if (configdof.n_vol() != shared_supercell->volume()) {
    std::stringstream msg;
    msg << "Error reading Configuration from JSON: "
        << "supercell size (" << shared_supercell->volume() << ") and "
        << "configdof size (" << configdof.n_vol() << ") are inconsistent.";
    msg << "supercell_name: " << json["supercell_name"].get<std::string>();
    parser.error.insert(msg.str());
  }
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  return notstd::make_unique<Configuration>(shared_supercell, configdof);
}

/// Insert Configuration to JSON
///
/// Format:
/// \code
/// {
///   "supercell_name": <string>,
///   "transformation_matrix_to_supercell": <3x3 array of integer>,
///   "dof": <JSON object representation of ConfigDoF>
/// }
/// \endcode
///
/// Note:
/// - "supercell_name": string
///       The name of the supercell. If the supercell is in canonical form, the
///       name has the format "SCELV_A_B_C_D_E_F", where V is the integer
///       supercell volume (multiple of the primitive cell), and "A_B_C_D_E_F"
///       are integer values of the hermite normal form of the
///       "transformation_matrix_to_supercell". If the supercell is not in
///       canonical form, the name has the format "SCELV_A_B_C_D_E_F.FG_INDEX",
///       where FG_INDEX is the index in the prim factor group of an operation
///       that transforms the canonical form to the particular supercell.
/// - "transformation_matrix_to_supercell": 3x3 array of integer
///       The integer transformation matrix T such that S = P*T, where P is a
///       3x3 column vector matrix of the prim lattice vectors and S is a 3x3
///       column vector matrix of the supercell lattice vectors.
/// - "dof": JSON object representation of ConfigDoF
///       See `to_json(ConfigDoF const &configdof, jsonParser &json)` for
///       details.
jsonParser &to_json(Configuration const &configuration, jsonParser &json) {
  if (!json.is_obj()) {
    throw std::runtime_error(
        "Error inserting configuration to json: not an object");
  }
  Supercell const &supercell = configuration.supercell();
  SupercellSymInfo const &sym_info = supercell.sym_info();
  json["supercell_name"] = supercell.name();
  json["transformation_matrix_to_supercell"] =
      sym_info.transformation_matrix_to_super();
  json["dof"] = configuration.configdof();
  // TODO: include properties, cache, source?

  return json;
}

/// Read Configuration from JSON
///
/// Note:
/// - See `to_json(Configuration const &configuration, jsonParser &json)` for
///   expected format.
/// - The attribute "transformation_matrix_to_super" must be present. If
///   "supercell_name" is present it is ignored.
/// - The constructed Configuration is not necessarily canonical, nor is its
///   Supercell.
/// - The constructed Configuration has a shared Supercell
///   (`std::shared_ptr<Supercell const>`) that is not owned by any Supercell
///   database.
/// - Use `DB::in_canonical_supercell` to make the canonical equivalent
///   configuration with the canonical supercell from the Supercell database,
///   without inserting the configuration in the Configuration database.
/// - Use `DB::make_canonical_and_insert` to make the equivalent configuration
///   with the canonical supercell from the Supercell database and insert the
///   canonical equivalent configuration in the Configuration database.
///
void from_json(Configuration &configuration, jsonParser const &json) {
  configuration = jsonConstructor<Configuration>::from_json(
      json, configuration.supercell().shared_prim());
  // TODO: include properties, cache, source?
}

}  // namespace CASM
