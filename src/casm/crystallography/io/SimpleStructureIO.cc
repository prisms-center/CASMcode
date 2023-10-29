#include "casm/crystallography/io/SimpleStructureIO.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/global/enum/json_io.hh"
#include "casm/global/errors.hh"

namespace {
using namespace CASM;

/// Read SimpleStructure from legacy `properties.calc.json` JSON format
///
/// Legacy `properties.calc.json` format, with prefix="relaxed":
/// \code
/// {
///   "atom_type\": [
///     "A",
///     "B"
///   ],
///   "atoms_per_type": [
///     1,
///     2
///   ],
///   "coordinate_mode": "direct",
///   "relaxed_basis": [
///     [0.6666667, 0.6666667, 0.6666667],
///     [0.00255632, 0.99488736, 0.00255632],
///     [0.33077698, 0.33844594, 0.33077698]
///   ],
///   "relaxed_energy\": -16.27773537,
///   "relaxed_forces\": [
///     [0.0, 0.0, 0.0],
///     [0.0, 0.00987362, -0.00987362],
///     [0.0, -0.00987362, 0.00987362]
///   ],
///   "relaxed_lattice\": [
///     [0.0, 1.9174843, 1.9174843],
///     [1.61158655, -1.88219884, 3.79968315],
///     [3.22317311, 0.0, 0.0]
///   ]
/// }
/// \endcode
///
/// Note:
/// - The properties "relaxed_mag_basis" and "relaxed_magmom" are no longer
/// supported. Now
///   magnetic properties are supported with the "magspin" flavors (see
///   AnisoValTraits documentation). If they are encountered they are silently
///   ignored.

static void _from_json_legacy(xtal::SimpleStructure &simple_structure,
                              const jsonParser &json) {
  try {
    COORD_TYPE coordinate_mode = json["coordinate_mode"].get<COORD_TYPE>();
    simple_structure.lat_column_mat =
        json["relaxed_lattice"].get<Eigen::Matrix3d>().transpose();

    // Input coordinate type to cartesian coordinate transformation matrix:
    //   cart_coord = coord_to_cart_mat * input_type_coord
    Eigen::Matrix3d to_cartesian_matrix = Eigen::Matrix3d::Identity();
    if (coordinate_mode == FRAC) {
      to_cartesian_matrix = simple_structure.lat_column_mat;
    }

    if (json.contains("relaxed_energy")) {
      simple_structure.properties["energy"] =
          json["relaxed_energy"].get<Eigen::MatrixXd>();
    }
    if (json.contains("relaxed_forces")) {
      simple_structure.atom_info.properties["force"] =
          json["relaxed_forces"].get<Eigen::MatrixXd>().transpose();
    }
    // note: silently ignoring legacy "relaxed_mag_basis" & "relaxed_magmom"

    // Read atom types
    std::vector<Index> atoms_per_type =
        json["atoms_per_type"].get<std::vector<Index> >();
    std::vector<std::string> atom_type =
        json["atom_type"].get<std::vector<std::string> >();
    for (Index i = 0; i < atoms_per_type.size(); ++i) {
      for (Index j = 0; j < atoms_per_type[i]; ++j) {
        simple_structure.atom_info.names.push_back(atom_type[i]);
      }
    }

    // Read basis coordinates
    simple_structure.atom_info.coords =
        to_cartesian_matrix *
        json["relaxed_basis"].get<Eigen::MatrixXd>().transpose();

  } catch (const std::exception &ex) {
    std::stringstream ss;
    ss << "Error parsing SimpleStructure from (legacy) JSON object. "
       << "One or more tags were improperly specified:\n"
       << ex.what();
    throw libcasm_runtime_error(ss.str());
  }
  return;
}

/// Output SimpleStructure atom_info and mol_info names to JSON
///
/// Note: "permute" enables skipping species that are excluded from output. It
/// is populated with indices of sites that should be included.
static jsonParser &_types_to_json(std::vector<std::string> const &names,
                                  jsonParser &json, std::string field_name,
                                  std::set<std::string> const &excluded_species,
                                  std::vector<Index> &permute) {
  jsonParser &names_json = json[field_name].put_array();
  for (Index i = 0; i < names.size(); ++i) {
    if (excluded_species.count(names[i])) continue;
    names_json.push_back(names[i]);
    permute.push_back(i);
  }
  return json;
}

/// Output SimpleStructure atom_info and mol_info coordinates to JSON
///
/// Note: "permute" enables skipping species that are excluded from output. It
/// includes only columns of the coords matrix that should be included.
static jsonParser &_coords_to_json(
    Eigen::MatrixXd const &coords, jsonParser &json, std::string field_name,
    std::vector<Index> const &permute,
    Eigen::Matrix3d const &to_coord_mode_matrix) {
  jsonParser &tjson = json[field_name].put_array();
  for (Index i : permute) {
    tjson.push_back(to_coord_mode_matrix * coords.col(i),
                    jsonParser::as_array());
  }
  return json;
}

/// Output SimpleStructure atom_info and mol_info properties to JSON
///
/// Note: "permute" enables skipping species that are excluded from output. It
/// includes only columns of the coords matrix that should be included.
static jsonParser &_properties_to_json(
    std::map<std::string, Eigen::MatrixXd> const &properties, jsonParser &json,
    std::string field_name, std::vector<Index> const &permute) {
  for (auto const &dof : properties) {
    jsonParser &tjson = json[field_name][dof.first]["value"].put_array();
    for (Index i : permute) {
      tjson.push_back(dof.second.col(i), jsonParser::as_array());
    }
  }
  return json;
}

/// Read SimpleStructure atom_info and mol_info type names from JSON
static void _types_from_json(std::vector<std::string> &names,
                             jsonParser const &json, std::string field_name) {
  if (json.contains(field_name)) {
    from_json(names, json[field_name]);
  }
}

/// Read SimpleStructure atom_info and mol_info coordinates from JSON
static void _coords_from_json(Eigen::MatrixXd &coords, jsonParser const &json,
                              std::string field_name,
                              Eigen::Matrix3d const &to_cartesian_matrix) {
  if (json.contains(field_name)) {
    coords = to_cartesian_matrix *
             json[field_name].get<Eigen::MatrixXd>().transpose();
  }
}

/// Read SimpleStructure global, atom_info, and mol_info properties from JSON
static void _properties_from_json(
    std::map<std::string, Eigen::MatrixXd> &properties, jsonParser const &json,
    std::set<std::string> allowed_field_names) {
  for (std::string const &field : allowed_field_names) {
    auto it = json.find(field);
    if (it != json.end()) {
      for (auto it2 = it->begin(); it2 != it->end(); ++it2) {
        Index min_data_size = 1;
        bool missing_data = false;
        if ( (*it2).size() < min_data_size ) { // no "value" field
          missing_data = true;
        }
        else if ( (*it2)["value"].size() < min_data_size ) { // empty "value"
          missing_data = true;
        }
        if (missing_data) {
          if (it2.name() == "force") {
            std::cout << "no force data (ignoring)\n";
            continue;
          }
          else {
            throw std::runtime_error(it2.name());
          }
        }
        properties[it2.name()] =
            (*it2)["value"].get<Eigen::MatrixXd>().transpose();
      }
    }
  }
}

/// Read SimpleStructure from JSON
///
/// \param simple_structure xtal::SimpleStructure to read from JSON
/// \param json A JSON object, from which the xtal::SimpleStructure JSON is
///     read. Expects JSON formatted as documented for `to_json` for
///     xtal::SimpleStructure.
///
static void _from_json_current(xtal::SimpleStructure &simple_structure,
                               const jsonParser &json) {
  COORD_TYPE coordinate_mode = json["coordinate_mode"].get<COORD_TYPE>();

  if (json.contains("lattice_vectors")) {
    simple_structure.lat_column_mat =
        json["lattice_vectors"].get<Eigen::Matrix3d>().transpose();
  } else {
    throw std::runtime_error(
        "Error reading xtal::SimpleStructure: \"lattice_vectors\" not found.");
  }

  // Input coordinate mode to cartesian coordinate transformation matrix:
  //   cart_coord = to_cartesian_matrix * input_mode_coord
  Eigen::Matrix3d to_cartesian_matrix;
  to_cartesian_matrix.setIdentity();
  if (coordinate_mode == FRAC) {
    to_cartesian_matrix = simple_structure.lat_column_mat;
  }

  // Read global properties (if exists)
  ::_properties_from_json(simple_structure.properties, json,
                          {"global_properties"});

  // Read atom_info (if exists)
  auto &atom_info = simple_structure.atom_info;
  ::_types_from_json(atom_info.names, json, "atom_type");
  ::_coords_from_json(atom_info.coords, json, "atom_coords",
                      to_cartesian_matrix);
  ::_properties_from_json(atom_info.properties, json, {"atom_properties"});

  // Read mol_info (if exists)
  auto &mol_info = simple_structure.mol_info;
  ::_types_from_json(mol_info.names, json, "mol_type");
  ::_coords_from_json(mol_info.coords, json, "mol_coords", to_cartesian_matrix);
  ::_properties_from_json(mol_info.properties, json, {"mol_properties"});
}

/// Check consistency of <atom/mol>_info.names and coords and properties size
///
/// Note: Will throw if mismatched sizes
static void _check_sizes(xtal::SimpleStructure::Info const &info,
                         std::string info_type) {
  if (info.names.size() != info.coords.cols()) {
    std::stringstream ss;
    ss << info_type << ".coords.cols(): " << info.coords.cols()
       << " != " << info_type << ".names.size(): " << info.names.size();
    throw libcasm_runtime_error(ss.str());
  }
  // Check consistency of <atom/mol>_info.names and properties size
  for (auto const &property_pair : info.properties) {
    if (info.names.size() != property_pair.second.cols()) {
      std::stringstream ss;
      ss << info_type << ".properties[" << property_pair.first
         << "\"].cols(): " << property_pair.second.cols() << " != " << info_type
         << ".names.size(): " << info.names.size();
      throw libcasm_runtime_error(ss.str());
    }
  }
}

}  // namespace

namespace CASM {

/// Output SimpleStructure to JSON
///
/// \param simple_structure xtal::SimpleStructure to output to JSON
/// \param json A JSON object, into which the xtal::SimpleStructure JSON will
/// output. This function
///        will overwrite existing attributes with the same name, but will not
///        erase any other existing attributes.
/// \param excluded_species Names of any molecular or atomic species that should
/// not be included
///        in the output
/// \param coordinate_mode COORD_TYPE (FRAC or CART) for output coordinates
///
/// Note:
/// - This is the expected format for `properties.calc.json`
///
/// Expected output:
/// \code
/// {
///   "coordinate_mode": ("Cartesian" (default), or "Direct" depending on
///   'mode') "atom_type": [  // array of atom type names
///     <atom type name>,
///     <atom type name>,
///     ...
///   ],
///   "mol_type": [ // array of molecule type names
///     <molecule type name>,
///     <molecule type name>,
///     ...
///   ],
///   "lattice_vectors": [
///      [<first lattice vector>],
///      [<second lattice vector>],
///      [<third lattice vector>]
///   ],
///   "global_properties": { // corresponds to simple_struc.properties
///      <property name>: {
///        "value": [<property vector>]
///      },
///      ...
///   }
///   "atom_properties": { // corresponds to simple_struc.atom_info.properties
///      <property name>: {
///        "value": [
///          [<property vector for site 0>],
///          [<property vector for site 1>],
///          ...
///        ]
///      },
///      ...
///   },
///   "mol_properties": { // corresponds to simple_struc.mol_info.properties
///      <property name>: {
///        "value": [
///          [<property vector for site 0>],
///          [<property vector for site 1>],
///          ...
///      },
///      ...
///   },
///   "atom_coords": [ // atom coordinates, according to "coordinate_mode"
///      [<coordinate for site 0>],
///      [<coordinate for site 1>],
///      ...
///   ],
///   "mol_coords": [ // molecule coordinates, according to "coordinate_mode"
///      [<coordinate for site 0>],
///      [<coordinate for site 1>],
///      ...
///   ],
///   ...
///   }
/// \endcode
///
jsonParser &to_json(xtal::SimpleStructure const &simple_structure,
                    jsonParser &json,
                    std::set<std::string> const &excluded_species,
                    COORD_TYPE coordinate_mode) {
  // Matrix to transform Cartesian coordinates to 'coordinate_mode' coordinates
  Eigen::Matrix3d to_coord_mode_matrix = Eigen::Matrix3d::Identity();
  if (coordinate_mode == FRAC) {
    json["coordinate_mode"] = FRAC;
    to_coord_mode_matrix = simple_structure.lat_column_mat.inverse();
  } else {
    json["coordinate_mode"] = CART;
  }

  json["lattice_vectors"] = simple_structure.lat_column_mat.transpose();

  // Output global properties
  for (auto const &dof : simple_structure.properties) {
    to_json_array(dof.second, json["global_properties"][dof.first]["value"]);
  }

  // Note: _types_to_json checks excluded_species and populates
  // "atom/mol_permute" with indices of
  //   sites that are included. This is used to only output the coordinates and
  //   properties of the correct sites in _properties_to_json and
  //   _coords_to_json

  // Output mol_info
  auto const &atom_info = simple_structure.atom_info;
  std::vector<Index> atom_permute;
  ::_types_to_json(atom_info.names, json, "atom_type", excluded_species,
                   atom_permute);
  ::_coords_to_json(atom_info.coords, json, "atom_coords", atom_permute,
                    to_coord_mode_matrix);
  ::_properties_to_json(atom_info.properties, json, "atom_properties",
                        atom_permute);

  // Output mol_info
  auto const &mol_info = simple_structure.mol_info;
  std::vector<Index> mol_permute;
  ::_types_to_json(mol_info.names, json, "mol_type", excluded_species,
                   mol_permute);
  ::_coords_to_json(mol_info.coords, json, "mol_coords", mol_permute,
                    to_coord_mode_matrix);
  ::_properties_to_json(mol_info.properties, json, "mol_properties",
                        mol_permute);

  return json;
}

//***************************************************************************

/// Read SimpleStructure from JSON
///
/// \param simple_structure xtal::SimpleStructure to read from JSON
/// \param json A JSON object, from which the xtal::SimpleStructure JSON is
/// read. Expects JSON
///        formatted either i) as documented for `to_json` for
///        xtal::SimpleStructure or ii) if the "atoms_per_type" attribute
///        exists, using the legacy `properties.calc.json` format as documented
///        for `from_json_legacy`.
///
/// Note:
/// - This is the expected format for `properties.calc.json`
///
void from_json(xtal::SimpleStructure &simple_structure,
               const jsonParser &json) {
  try {
    // For backwards read compatibility:
    if (json.contains("atoms_per_type")) {
      ::_from_json_legacy(simple_structure, json);
    } else {
      ::_from_json_current(simple_structure, json);
    }

    // Check consistency of <atom/mol>_info.names and coords and properties size
    ::_check_sizes(simple_structure.atom_info, "atom_info");
    ::_check_sizes(simple_structure.mol_info, "mol_info");

  } catch (const std::exception &ex) {
    std::stringstream ss;
    ss << "Error parsing SimpleStructure from JSON object. "
       << "One or more tags were improperly specified:\n"
       << ex.what();
    throw std::runtime_error(ss.str());
  }
}
}  // namespace CASM
