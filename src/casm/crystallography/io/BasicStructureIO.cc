#include "casm/crystallography/io/BasicStructureIO.hh"

#include <boost/filesystem.hpp>

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/SafeOfstream.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/AnisoValTraits.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/crystallography/io/DoFSetIO.hh"
#include "casm/global/enum/json_io.hh"
#include "casm/misc/ParsingDictionary.hh"

namespace CASM {

jsonParser const &from_json(xtal::SpeciesProperty &_attr,
                            jsonParser const &json) {
  _attr.set_value(json["value"].get<Eigen::VectorXd>());
  return json;
}

//****************************************************

jsonParser &to_json(xtal::SpeciesProperty const &_attr, jsonParser &json) {
  json.put_obj();
  to_json_array(_attr.value(), json["value"]);
  return json;
}

//****************************************************
jsonParser &to_json(xtal::AtomPosition const &apos, jsonParser &json,
                    Eigen::Ref<const Eigen::Matrix3d> const &c2f_mat) {
  json.put_obj();
  to_json_array(c2f_mat * apos.cart(), json["coordinate"]);
  json["name"] = apos.name();
  if (apos.properties().size()) json["properties"] = apos.properties();
  return json;
}

//****************************************************

void from_json(xtal::AtomPosition &apos, const jsonParser &json,
               Eigen::Ref<const Eigen::Matrix3d> const &f2c_mat,
               ParsingDictionary<AnisoValTraits> const &_aniso_val_dict) {
  apos = json.get<xtal::AtomPosition>(f2c_mat, _aniso_val_dict);

  return;
}

//****************************************************

xtal::AtomPosition jsonConstructor<xtal::AtomPosition>::from_json(
    const jsonParser &json, Eigen::Matrix3d const &f2c_mat,
    ParsingDictionary<AnisoValTraits> const &_aniso_val_dict) {
  std::string _name;
  Eigen::Vector3d _pos(0., 0., 0.);
  std::map<std::string, xtal::SpeciesProperty> attr_map;
  if (json.is_obj()) {
    _name = json["name"].get<std::string>();
    if (json.contains("coordinate")) {
      _pos = f2c_mat * json["coordinate"].get<Eigen::Vector3d>();
      // std::cout << "f2c_mat: \n" << f2c_mat << "\n";
      // std::cout << "_pos: " << _pos.transpose() << "\n";
    }
    if (json.contains("properties")) {
      auto it = json["properties"].cbegin(), end_it = json["properties"].cend();
      for (; it != end_it; ++it) {
        auto result_pair =
            attr_map.emplace(it.name(), _aniso_val_dict.lookup(it.name()));
        CASM::from_json(result_pair.first->second, *it);
      }
    }

  } else if (json.is_string()) {
    _name = json.get<std::string>();
  } else
    throw std::runtime_error(
        "Invalid JSON input encountered. Unable to parse AtomPosition "
        "object.\n");

  xtal::AtomPosition result(_pos, _name);
  result.set_properties(attr_map);
  return result;
}

//****************************************************
//   Write Molecule to json.
//****************************************************

jsonParser &to_json(xtal::Molecule const &mol, jsonParser &json,
                    Eigen::Ref<const Eigen::Matrix3d> const &c2f_mat) {
  json.put_obj();
  CASM::to_json(mol.atoms(), json["atoms"], c2f_mat);
  json["name"] = mol.name();
  if (mol.properties().size()) json["properties"] = mol.properties();

  return json;
}

//****************************************************
//
//    Read Molecule from json.
//
//****************************************************

void from_json(xtal::Molecule &mol, const jsonParser &json,
               Eigen::Ref<const Eigen::Matrix3d> const &f2c_mat,
               ParsingDictionary<AnisoValTraits> const &_aniso_val_dict) {
  if (json.contains("atoms")) {
    std::vector<xtal::AtomPosition> _atoms;
    CASM::from_json(_atoms, json["atoms"], f2c_mat, _aniso_val_dict);
    mol.set_atoms(_atoms);
  }

  std::map<std::string, xtal::SpeciesProperty> attr_map;
  if (json.contains("properties")) {
    auto it = json["properties"].cbegin(), end_it = json["properties"].cend();
    for (; it != end_it; ++it) {
      auto result_pair =
          attr_map.emplace(it.name(), _aniso_val_dict.lookup(it.name()));
      from_json(result_pair.first->second, *it);
    }
  }
  mol.set_properties(attr_map);

  // jsonParser tjson;
  // to_json(mol,tjson,f2c_mat.inverse());
  // std::cout << "Read Molecule :\n"<< json << "\n";
}

//****************************************************
//
//****************************************************

xtal::Molecule jsonConstructor<xtal::Molecule>::from_json(
    const jsonParser &json, Eigen::Ref<const Eigen::Matrix3d> const &f2c_mat,
    ParsingDictionary<AnisoValTraits> const &_aniso_val_dict) {
  return json.get<xtal::Molecule>(f2c_mat, _aniso_val_dict);
}

//****************************************************

xtal::Site jsonConstructor<xtal::Site>::from_json(
    const jsonParser &json, xtal::Lattice const &_home, COORD_TYPE coordtype,
    std::map<std::string, xtal::Molecule> const &mol_map,
    std::vector<std::vector<std::string>> &unique_names,
    ParsingDictionary<AnisoValTraits> const &_aniso_val_dict) {
  xtal::Site result(_home);
  std::vector<std::string> site_unique_names;
  CASM::from_json(result, json, _home, coordtype, mol_map, site_unique_names,
                  _aniso_val_dict);
  unique_names.push_back(site_unique_names);
  return result;
}

//****************************************************

jsonParser &to_json(const xtal::Site &site, jsonParser &json,
                    COORD_TYPE coordtype) {
  json.put_obj();

  // class xtal::Site : public Coordinate
  if (coordtype == FRAC)
    to_json_array(site.frac(), json["coordinate"]);
  else
    to_json_array(site.cart(), json["coordinate"]);

  // Occupant_DoF<Molecule> occupant_dof;
  Eigen::Matrix3d c2f = Eigen::Matrix3d::Identity();
  // change this to use FormatFlag
  if (coordtype == FRAC) c2f = site.home().inv_lat_column_mat();
  CASM::to_json(site.occupant_dof(), json["occupants"], c2f);

  if (site.dofs().size()) {
    json["dofs"] = site.dofs();
  }

  // Index m_label
  if (valid_index(site.label())) json["label"] = site.label();

  return json;
}

//****************************************************

void from_json(xtal::Site &site, const jsonParser &json,
               xtal::Lattice const &_home, COORD_TYPE coordtype,
               std::map<std::string, xtal::Molecule> const &mol_map,
               std::vector<std::string> &site_unique_names,
               ParsingDictionary<AnisoValTraits> const &_aniso_val_dict) {
  site.set_lattice(_home, coordtype);
  if (coordtype == FRAC)
    site.frac() = json["coordinate"].get<Eigen::Vector3d>();
  else
    site.cart() = json["coordinate"].get<Eigen::Vector3d>();

  // Index m_label -- must be greater than zero
  Index _label = -1;
  if (json.contains("label")) {
    CASM::from_json(_label, json["label"]);
    if (!valid_index(_label))
      throw std::runtime_error("JSON specification of site has {\"label\" : " +
                               std::to_string(_label) +
                               "}, but \"label\" must be greater than 0.\n");
  }
  site.set_label(_label);

  // Local continuous dofs

  std::map<std::string, xtal::SiteDoFSet> _dof_map;
  if (json.contains("dofs")) {
    auto it = json["dofs"].begin(), end_it = json["dofs"].end();
    for (; it != end_it; ++it) {
      if (_dof_map.count(it.name()))
        throw std::runtime_error(
            "Error parsing local field \"dofs\" from JSON. DoF type " +
            it.name() + " cannot be repeated.");

      try {
        /* _dof_map.emplace(std::make_pair(it.name(),
         * it->get<xtal::DoFSet>(_aniso_val_dict.lookup(it.name()))));
         */
        _dof_map.emplace(std::make_pair(
            it.name(),
            it->get<xtal::SiteDoFSet>(_aniso_val_dict.lookup(it.name()))));
        /* _dof_map.emplace(std::make_pair(it.name(),
         * from_json<xtal::SiteDoFSet>(*it))); */
      } catch (std::exception &e) {
        throw std::runtime_error(
            "Error parsing local field \"dofs\" from JSON. Failure for DoF "
            "type " +
            it.name() + ": " + e.what());
      }
    }
  }
  site.set_dofs(_dof_map);

  std::vector<xtal::Molecule> t_occ;
  std::string occ_key;
  if (json.contains("occupants"))
    occ_key = "occupants";
  else if (json.contains("occupant_dof"))
    occ_key = "occupant_dof";

  if (!occ_key.empty()) {
    for (std::string const &occ :
         json[occ_key].get<std::vector<std::string>>()) {
      // std::cout << "CREATING OCCUPANT " << occ << "\n";
      // Have convenience options for properties like magnetic moment, etc?
      site_unique_names.push_back(occ);
      auto it = mol_map.find(occ);
      if (it != mol_map.end())
        t_occ.push_back(it->second);
      else
        t_occ.push_back(xtal::Molecule::make_atom(occ));
    }
  }
  // std::cout << "t_occ.size() = " << t_occ.size() << "\n";
  if (t_occ.empty()) t_occ.push_back(xtal::Molecule::make_unknown());
  site.set_allowed_occupants(t_occ);
}

/// Read BasicStructure, detecting supported file types
///
/// \param filename Path to BasicStructure input file
/// \param xtal_tol Lattice tolerance to use
/// \param _aniso_val_dict Can be used for reading custom AnisoValTraits
xtal::BasicStructure read_prim(
    fs::path filename, double xtal_tol,
    ParsingDictionary<AnisoValTraits> const *_aniso_val_dict) {
  std::string prim_file_type;
  return read_prim(filename, xtal_tol, _aniso_val_dict, prim_file_type);
}

/// Read BasicStructure, detecting supported file types
///
/// \param filename Path to BasicStructure input file
/// \param xtal_tol Lattice tolerance to use
/// \param _aniso_val_dict Can be used for reading custom AnisoValTraits
/// \param prim_file_type When successful, set with input file type. One of
///     "vasp", "json".
xtal::BasicStructure read_prim(
    fs::path filename, double xtal_tol,
    ParsingDictionary<AnisoValTraits> const *_aniso_val_dict,
    std::string &prim_file_type) {
  jsonParser json;
  filename = fs::absolute(filename);
  if (!fs::exists(filename)) {
    throw std::invalid_argument("Error reading prim from file '" +
                                filename.string() + "'. Does not exist.");
  }
  std::ifstream f(filename.string());
  if (!f) {
    throw std::invalid_argument("Error reading prim from file '" +
                                filename.string() +
                                "'. Does not contain valid input.");
  }
  while (!f.eof() && std::isspace(f.peek())) {
    f.get();
  }

  // TODO: Use functions. One for checking file type, one for reading json, one
  // for reading vasp.
  // Check if JSON
  if (f.peek() == '{') {
    try {
      prim_file_type = "json";
      json = jsonParser(f);
    } catch (std::exception const &ex) {
      std::stringstream err_msg;
      err_msg << "Error reading prim from JSON file '" << filename
              << "':" << std::endl
              << ex.what();
      throw std::invalid_argument(err_msg.str());
    }
    return read_prim(json, xtal_tol, _aniso_val_dict);
  }
  // else tread as vasp-like file
  xtal::BasicStructure prim;
  try {
    prim_file_type = "vasp";
    prim = xtal::BasicStructure::from_poscar_stream(f);
  } catch (std::exception const &ex) {
    std::stringstream err_msg;
    err_msg << "Error reading prim from VASP-formatted file '" << filename
            << "':" << std::endl
            << ex.what();
    throw std::invalid_argument(err_msg.str());
  }
  return prim;
}

/// \brief Read prim.json
xtal::BasicStructure read_prim(
    const jsonParser &json, double xtal_tol,
    ParsingDictionary<AnisoValTraits> const *_aniso_val_dict) {
  ParsingDictionary<AnisoValTraits> default_aniso_val_dict =
      make_parsing_dictionary<AnisoValTraits>();
  if (_aniso_val_dict == nullptr) _aniso_val_dict = &default_aniso_val_dict;

  // read lattice
  Eigen::Matrix3d latvec_transpose;

  try {
    from_json(latvec_transpose, json["lattice_vectors"]);
  } catch (std::exception &e) {
    log() << e.what() << std::endl;
    throw std::runtime_error(
        "Error parsing global field \"lattice_vectors\" from prim JSON.");
  }

  xtal::Lattice lat(latvec_transpose.transpose(), xtal_tol);

  // create prim using lat
  xtal::BasicStructure prim(lat);

  // read title
  try {
    prim.set_title(json["title"].get<std::string>());
  } catch (std::exception &e) {
    log() << e.what() << std::endl;
    throw std::runtime_error(
        "Error parsing global field \"title\" from prim JSON.");
  }

  Eigen::Vector3d vec;

  // Global DoFs
  try {
    std::map<std::string, xtal::DoFSet> _dof_map;
    if (json.contains("dofs")) {
      auto it = json["dofs"].begin(), end_it = json["dofs"].end();
      for (; it != end_it; ++it) {
        if (_dof_map.count(it.name()))
          throw std::runtime_error(
              "Error parsing global field \"dofs\" from prim JSON. DoF type " +
              it.name() + " cannot be repeated.");

        try {
          // TODO: Am I messing something up here? Why was it constructed so
          // weird before?
          /* _dof_map.emplace(std::make_pair(it.name(),
           * it->get<xtal::DoFSet>(_aniso_val_dict.lookup(it.name()))));
           */
          _dof_map.emplace(std::make_pair(
              it.name(),
              it->get<xtal::DoFSet>(_aniso_val_dict->lookup(it.name()))));
        } catch (std::exception &e) {
          log() << e.what() << std::endl;
          throw std::runtime_error(
              "Error parsing global field \"dofs\" from prim JSON. Failure for "
              "DoF type " +
              it.name() + ": " + e.what());
        }
      }
      prim.set_global_dofs(_dof_map);
    }
  } catch (std::exception &e) {
    log() << e.what() << std::endl;
    throw std::runtime_error(
        "Error parsing global field \"dofs\" from prim JSON.");
  }

  // read basis coordinate mode
  COORD_TYPE mode;
  try {
    from_json(mode, json["coordinate_mode"]);
  } catch (std::exception &e) {
    log() << e.what() << std::endl;
    throw std::runtime_error(
        "Error parsing global field \"coordinate_mode\" from prim JSON.");
  }

  // Molecules
  std::map<std::string, xtal::Molecule> mol_map;
  Eigen::Matrix3d f2c;
  if (mode == FRAC)
    f2c = lat.lat_column_mat();
  else
    f2c.setIdentity();

  try {
    if (json.contains("species")) {
      auto it = json["species"].begin();
      auto it_end = json["species"].end();
      for (; it != it_end; ++it) {
        std::string chem_name = it.name();
        it->get_if(chem_name, "name");
        // log() << "chem_name: " << chem_name << "\n";
        auto mol_it =
            mol_map.emplace(it.name(), xtal::Molecule(chem_name)).first;
        from_json(mol_it->second, *it, f2c, *_aniso_val_dict);
      }
    }
  } catch (std::exception &e) {
    log() << e.what() << std::endl;
    throw std::runtime_error(
        "Error parsing global field \"species\" from prim JSON.");
  }

  try {
    // read basis sites
    std::vector<std::vector<std::string>> unique_names;
    for (jsonParser const &bjson : json["basis"])
      prim.push_back(bjson.get<xtal::Site>(prim.lattice(), mode, mol_map,
                                           unique_names, *_aniso_val_dict));
    prim.set_unique_names(unique_names);
  } catch (std::exception &e) {
    log() << e.what() << std::endl;
    throw std::runtime_error(
        "Error parsing global field \"basis\" from prim JSON.");
  }
  return prim;
}

/// \brief Write prim.json to file
void write_prim(const xtal::BasicStructure &prim, fs::path filename,
                COORD_TYPE mode, bool include_va) {
  SafeOfstream outfile;
  outfile.open(filename);

  jsonParser json;
  write_prim(prim, json, mode, include_va);
  json.print(outfile.ofstream());

  outfile.close();
}

/// \brief Write prim.json as JSON
void write_prim(const xtal::BasicStructure &prim, jsonParser &json,
                COORD_TYPE mode, bool include_va) {
  // TODO: Use functions. One for each type that's getting serialized.

  json = jsonParser::object();

  json["title"] = prim.title();

  json["lattice_vectors"] = prim.lattice().lat_column_mat().transpose();

  if (mode == COORD_DEFAULT) {
    mode = xtal::COORD_MODE::CHECK();
  }

  Eigen::Matrix3d c2f_mat;
  if (mode == FRAC) {
    c2f_mat = prim.lattice().inv_lat_column_mat();
    json["coordinate_mode"] = FRAC;
  } else if (mode == CART) {
    c2f_mat.setIdentity();
    json["coordinate_mode"] = CART;
  }

  // Global DoFs
  for (auto const &_dof : prim.global_dofs()) {
    json["dofs"][_dof.first] = _dof.second;
  }

  std::vector<std::vector<std::string>> mol_names = prim.unique_names();
  if (mol_names.empty()) {
    mol_names = allowed_molecule_unique_names(prim);
  }
  jsonParser &bjson = (json["basis"].put_array());
  for (int i = 0; i < prim.basis().size(); i++) {
    xtal::Site site = prim.basis()[i];
    if (!include_va && site.occupant_dof().size() == 1 &&
        site.occupant_dof()[0].is_vacancy())
      continue;
    bjson.push_back(jsonParser::object());
    jsonParser &sjson = bjson[bjson.size() - 1];
    if (valid_index(site.label())) sjson["label"] = site.label();
    jsonParser &cjson = sjson["coordinate"].put_array();
    if (mode == FRAC) {
      cjson.push_back(site.frac(0));
      cjson.push_back(site.frac(1));
      cjson.push_back(site.frac(2));
    } else if (mode == CART) {
      cjson.push_back(site.cart(0));
      cjson.push_back(site.cart(1));
      cjson.push_back(site.cart(2));
    }

    if (site.dofs().size()) {
      sjson["dofs"] = site.dofs();
    }

    jsonParser &ojson =
        (sjson["occupants"] = jsonParser::array(site.occupant_dof().size()));

    for (int j = 0; j < site.occupant_dof().size(); j++) {
      ojson[j] = mol_names[i][j];
      if (site.occupant_dof()[j].name() != mol_names[i][j] ||
          !site.occupant_dof()[j].is_atomic()) {
        to_json(site.occupant_dof()[j], json["species"][mol_names[i][j]],
                c2f_mat);
      }
    }
  }
}

void from_json(xtal::BasicStructure &prim, jsonParser const &json,
               double xtal_tol,
               ParsingDictionary<AnisoValTraits> const *_aniso_val_dict) {
  prim = read_prim(json, xtal_tol, _aniso_val_dict);
}

jsonParser &to_json(const xtal::BasicStructure &prim, jsonParser &json,
                    COORD_TYPE mode, bool include_va) {
  write_prim(prim, json, mode, include_va);
  return json;
}

}  // namespace CASM
