#include "casm/clex/io/json/ConfigDoF_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

jsonParser &to_json(LocalContinuousConfigDoFValues const &_values,
                    jsonParser &_json) {
  to_json(_values.standard_values().transpose(), _json["values"]);
  return _json;
}

void from_json(LocalContinuousConfigDoFValues &_values,
               jsonParser const &_json) {
  Eigen::MatrixXd tval = _json["values"].get<Eigen::MatrixXd>().transpose();
  _values.from_standard_values(tval);
}

jsonParser &to_json(LocalDiscreteConfigDoFValues const &_values,
                    jsonParser &_json) {
  to_json_array(_values.values(), _json);
  return _json;
}

void from_json(LocalDiscreteConfigDoFValues &_values, jsonParser const &_json) {
  _values.resize_vol(_json.size() / _values.n_sublat());
  _values.values() = _json.get<LocalDiscreteConfigDoFValues::ValueType>();
}

jsonParser &to_json(GlobalContinuousConfigDoFValues const &_values,
                    jsonParser &_json) {
  to_json_array(_values.standard_values(), _json["values"]);
  return _json;
}

void from_json(GlobalContinuousConfigDoFValues &_values,
               jsonParser const &_json) {
  _values.from_standard_values(_json["values"].get<Eigen::VectorXd>());
}

std::unique_ptr<ConfigDoF> jsonMake<ConfigDoF>::make_from_json(
    const jsonParser &json, Structure const &prim) {
  auto result_ptr = notstd::make_unique<ConfigDoF>(
      (Index)prim.basis().size(), 0, global_dof_info(prim),
      local_dof_info(prim), prim.occupant_symrep_IDs(), prim.lattice().tol());
  result_ptr->from_json(json);
  return result_ptr;
}

ConfigDoF jsonConstructor<ConfigDoF>::from_json(const jsonParser &json,
                                                Structure const &prim) {
  ConfigDoF result{(Index)prim.basis().size(), 0,
                   global_dof_info(prim),      local_dof_info(prim),
                   prim.occupant_symrep_IDs(), prim.lattice().tol()};
  result.from_json(json);
  return result;
}

/// Insert ConfigDoF to JSON
///
/// Format:
/// \code
/// {
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
  if (configdof.occupation().size())
    to_json_array(configdof.occupation(), json["occ"]);
  if (!configdof.local_dofs().empty()) {
    json["local_dofs"] = configdof.local_dofs();
  }
  if (!configdof.global_dofs().empty()) {
    json["global_dofs"] = configdof.global_dofs();
  }

  return json;
}

/// Read ConfigDoF from JSON
void from_json(ConfigDoF &configdof, const jsonParser &json) {
  if (json.contains("occupation")) {
    // //For Backwards compatibility
    configdof.set_occupation(
        json["occupation"].get<LocalDiscreteConfigDoFValues>().values());
  } else if (json.contains("occ")) {
    configdof.set_occupation(
        json["occ"].get<LocalDiscreteConfigDoFValues>().values());
  } else {
    throw std::runtime_error(
        "JSON serialization of ConfigDoF must contain field \"occ\"\n");
  }

  if (json.contains("local_dofs")) {
    auto end_it = json["local_dofs"].end();
    for (auto it = json["local_dofs"].begin(); it != end_it; ++it) {
      CASM::from_json(configdof.local_dof(it.name()), *it);
    }
  }

  if (json.contains("global_dofs")) {
    auto end_it = json["global_dofs"].end();
    for (auto it = json["global_dofs"].begin(); it != end_it; ++it) {
      CASM::from_json(configdof.global_dof(it.name()), *it);
    }
  }
}

}  // namespace CASM
