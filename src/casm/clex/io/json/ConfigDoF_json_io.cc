#include "casm/clex/io/json/ConfigDoF_json_io.hh"

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/ConfigDoFTools.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

std::unique_ptr<ConfigDoF> jsonMake<ConfigDoF>::make_from_json(
    const jsonParser &json, Structure const &prim, Index volume) {
  auto result_ptr = make_unique_configdof(prim, volume);
  from_json(*result_ptr, json);
  return result_ptr;
}

ConfigDoF jsonConstructor<ConfigDoF>::from_json(const jsonParser &json,
                                                Structure const &prim,
                                                Index volume) {
  ConfigDoF result = make_configdof(prim, volume);
  CASM::from_json(result, json);
  return result;
}

/// Insert ConfigDoF to JSON
///
/// Format:
/// \code
/// {
///   "volume": <int>,
///   "occ": <array of integer>,
///   "global_dofs": {
///     <dof_name>: {
///       "values": <1d array of numbers, DoF values in standard basis>
///     },
///     ...
///   },
///   "local_dofs": {
///     <dof_name>: {
///       "values": <2d array of numbers, DoF values in standard basis>
///     },
///     ...
///   }
/// }
/// \endcode
///
/// Reminder about standard DoF basis vs prim DoF basis:
/// - Each type of DoF refers to an AnisoValTraits object which contains the
/// "standard" DoF basis
/// - The prim BasicStructure holds xtal::SiteDoFSet (representing local DoF)
/// and xtal::DoFSet
///   (representing global DoF) which contain an AnisoValTraits object and a
///   separate "prim" DoF basis containing a set of named basis vectors which
///   are denoted relative to the standard basis, allowing the user to specify
///   the DoFSet components, name them, and restrict DoF values to a particular
///   subspace.
/// - Examples of standard DoF basis specified by AnisoValTraits:
///   - "disp" -> (dx, dy, dz) -> displacement components relative to fixed
///   laboratory frame
///   - "strain" -> (e_xx, e_yy, e_zz, sqrt(2)*e_yz, sqrt(2)*e_xz, sqrt(2)*e_xy)
///   -> tensor elements
/// - ConfigDoF objects in memory store DoF values as coordinates in the prim
/// DoF basis
/// - ConfigDoF JSON representation stores DoF values as coordinates in the
/// standard DoF basis
///
/// Note:
/// - "occ":
///       The integer value for each site corresponding to which Molecule in the
///       Site::occupant_dof vector is occupying that site.
///
///       Example: supercell volume=3, prim basis size=2,
///
///           "occ": [
///               occ[1], // "occ" value on site sublattice 0, unit cell 0
///               occ[2], // "occ" value on site sublattice 0, unit cell 1
///               occ[3], // "occ" value on site sublattice 0, unit cell 2
///               occ[4], // "occ" value on site sublattice 1, unit cell 0
///               occ[5], // "occ" value on site sublattice 1, unit cell 1
///               occ[6]  // "occ" value on site sublattice 1, unit cell 2
///             ]
///
/// - "local_dofs":
///       The local DoF values are represented as a matrix, with each row
///       representing a site DoF value and each colum representing a component
///       of the DoF value:
///           number of cols = DoF dimension (i.e. 3 for "disp")
///           number of rows = (Supercell volume as multiple of the prim) *
///           (prim basis size).
///
///       Example: "disp" values, supercell volume=3, prim basis size=2
///
///           "local_dofs": {
///             "disp": {
///               "values": [
///                 [dx[1], dy[1], dz[1]], // "disp" values on site: sublattice
///                 0, unit cell 0 [dx[2], dy[2], dz[2]], // "disp" values on
///                 site: sublattice 0, unit cell 1 [dx[3], dy[3], dz[3]], //
///                 "disp" values on site: sublattice 0, unit cell 2 [dx[4],
///                 dy[4], dz[4]], // "disp" values on site: sublattice 1, unit
///                 cell 0 [dx[5], dy[5], dz[5]], // "disp" values on site:
///                 sublattice 1, unit cell 1 [dx[6], dy[6], dz[6]], // "disp"
///                 values on site: sublattice 1, unit cell 2
///               ]
///             }
///           }
///
/// - "global_dofs":
///       The global DoF values are represented as a vector of size equal to the
///       dimension of the DoF (i.e. 6 for "GLstrain").
///
///       Example: "GLstrain" values (any supercell volume and prim basis size)
///
///           "global_dofs": {
///             "GLstrain": {
///               "values": [Exx, Eyy, Ezz, sqrt(2)Exz, sqrt(2)Eyz, sqrt(2)Exy]
///             }
///           }
///
jsonParser &to_json(const ConfigDoF &configdof, jsonParser &json) {
  json = jsonParser::object();
  if (configdof.occupation().size()) {
    to_json_array(configdof.occupation(), json["occ"]);
  }
  if (!configdof.local_dofs().empty()) {
    for (auto const &local_dof : configdof.local_dofs()) {
      to_json(local_dof.second.standard_values().transpose(),
              json["local_dofs"][local_dof.first]["values"]);
    }
  }
  if (!configdof.global_dofs().empty()) {
    for (auto const &global_dof : configdof.global_dofs()) {
      to_json_array(global_dof.second.standard_values(),
                    json["global_dofs"][global_dof.first]["values"]);
    }
  }

  return json;
}

/// Read ConfigDoF from JSON
void from_json(ConfigDoF &configdof, const jsonParser &json) {
  auto &log = CASM::log();
  ParentInputParser parser{json};
  std::runtime_error error_if_invalid{
      "Error reading ConfigDoF from JSON input"};

  Eigen::VectorXi occ;
  // For Backwards compatibility
  if (parser.self.contains("occupation")) {
    parser.require(occ, "occupation");
  } else {
    parser.require(occ, "occ");
    try {
      configdof.set_occupation(occ);
    } catch (std::exception const &e) {
      parser.insert_error("occ", e.what());
    }
  }

  for (auto const &local_dof : configdof.local_dofs()) {
    Eigen::MatrixXd tvalues;
    fs::path values_path = fs::path{"local_dofs"} / local_dof.first / "values";
    parser.require(tvalues, values_path);
    try {
      configdof.local_dof(local_dof.first)
          .from_standard_values(tvalues.transpose());
    } catch (std::exception const &e) {
      parser.insert_error(values_path, e.what());
    }
  }

  for (auto const &global_dof : configdof.global_dofs()) {
    Eigen::VectorXd tvalues;
    fs::path values_path =
        fs::path{"global_dofs"} / global_dof.first / "values";
    parser.require(tvalues, values_path);
    try {
      configdof.global_dof(global_dof.first).from_standard_values(tvalues);
    } catch (std::exception const &e) {
      parser.insert_error(values_path, e.what());
    }
  }

  report_and_throw_if_invalid(parser, log, error_if_invalid);
}

}  // namespace CASM
