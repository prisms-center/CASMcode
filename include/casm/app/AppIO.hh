#ifndef CASM_AppIO
#define CASM_AppIO

#include "casm/crystallography/BasicStructure.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/clex/ChemicalReference.hh"

#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/json_io/clex.hh"

namespace CASM {

  // --- These functions are for the casm executable file I/O -----------
  //
  //  They should be moved into the app code eventually, but are currently needed here


  // --------- PrimIO Declarations --------------------------------------------------

  BasicStructure<Site> read_prim(fs::path filename);

  BasicStructure<Site> read_prim(const jsonParser &json);

  /// \brief Write prim.json to file
  void write_prim(const BasicStructure<Site> &prim, fs::path filename, COORD_TYPE mode);

  /// \brief Write prim.json as JSON
  void write_prim(const BasicStructure<Site> &prim, jsonParser &json, COORD_TYPE mode);


  // --------- PrimIO Definitions --------------------------------------------------

  inline BasicStructure<Site> read_prim(fs::path filename) {

    try {
      jsonParser json(filename);
      return read_prim(json);
    }
    catch(...) {
      std::cerr << "Error reading prim from " << filename << std::endl;
      /// re-throw exceptions
      throw;
    }
  }

  /// \brief Read prim.json
  inline BasicStructure<Site> read_prim(const jsonParser &json) {

    try {

      // read lattice
      Vector3<double> vec0, vec1, vec2;

      from_json(vec0, json["lattice_vectors"][0]);
      from_json(vec1, json["lattice_vectors"][1]);
      from_json(vec2, json["lattice_vectors"][2]);

      Lattice lat(vec0, vec1, vec2);

      // create prim using lat
      BasicStructure<Site> prim(lat);

      // read title
      from_json(prim.title, json["title"]);

      Vector3<double> vec;

      // read basis coordinate mode
      std::string coordinate_mode;
      from_json(coordinate_mode, json["coordinate_mode"]);
      COORD_TYPE mode;
      if(coordinate_mode == "Fractional") {
        mode = FRAC;
      }
      else if(coordinate_mode == "Direct") {
        mode = FRAC;
      }
      else if(coordinate_mode == "Cartesian") {
        mode = CART;
      }
      else {
        throw std::runtime_error(
          std::string(" Invalid: \"coordinate_mode\"\n") +
          "   Expected one of \"Fractional\", \"Direct\", or \"Cartesian\"");
      }

      // read basis sites
      Vector3<double> coord;

      for(int i = 0; i < json["basis"].size(); i++) {

        // read coordinate
        from_json(coord, json["basis"][i]["coordinate"]);
        Site site(prim.lattice());
        site(mode) = coord;

        // read atom occupant names
        Array<std::string> occ_name;
        from_json(occ_name, json["basis"][i]["occupant_dof"]);

        // fill site.site_occupant
        Array<Molecule> tocc;
        for(int i = 0; i < occ_name.size(); i++) {
          Molecule tMol(prim.lattice());
          tMol.name = occ_name[i];
          tMol.push_back(AtomPosition(0, 0, 0, occ_name[i], prim.lattice()));
          tocc.push_back(tMol);
        }
        site.set_site_occupant(MoleculeOccupant(tocc));
        site.set_occ_value(0);

        // add site to prim
        prim.basis.push_back(site);
      }

      return prim;
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }

  }

  /// \brief Write prim.json to file
  inline void write_prim(const BasicStructure<Site> &prim, fs::path filename, COORD_TYPE mode) {

    SafeOfstream outfile;
    outfile.open(filename);

    jsonParser json;
    write_prim(prim, json, mode);
    json.print(outfile.ofstream());

    outfile.close();
  }

  /// \brief Write prim.json as JSON
  inline void write_prim(const BasicStructure<Site> &prim, jsonParser &json, COORD_TYPE mode) {

    json = jsonParser::object();

    json["title"] = prim.title;

    json["lattice_vectors"] = jsonParser::array();
    json["lattice_vectors"].push_back(prim.lattice()[0]);
    json["lattice_vectors"].push_back(prim.lattice()[1]);
    json["lattice_vectors"].push_back(prim.lattice()[2]);

    if(mode == COORD_DEFAULT) {
      mode = COORD_MODE::CHECK();
    }

    if(mode == FRAC) {
      json["coordinate_mode"] = "Fractional";
    }
    else if(mode == CART) {
      json["coordinate_mode"] = "Cartesian";
    }

    json["basis"] = jsonParser::array(prim.basis.size());
    for(int i = 0; i < prim.basis.size(); i++) {
      json["basis"][i] = jsonParser::object();
      json["basis"][i]["coordinate"] = prim.basis[i](mode);
      json["basis"][i]["occupant_dof"] = jsonParser::array(prim.basis[i].site_occupant().size());

      for(int j = 0; j < prim.basis[i].site_occupant().size(); j++) {
        json["basis"][i]["occupant_dof"][j] = prim.basis[i].site_occupant()[j].name;
      }

    }
  }



  // --------- SymmetryIO Declarations --------------------------------------------------

  inline void write_symop(const SymOp &op, jsonParser &json, int cclass, int inv) {
    json = jsonParser::object();

    json["matrix"]["CART"] = op.get_matrix(CART);
    json["matrix"]["FRAC"] = op.get_matrix(FRAC);
    json["tau"]["CART"] = op.tau(CART);
    json["tau"]["FRAC"] = op.tau(FRAC);
    json["conjugacy_class"] = cclass;
    json["inverse"] = inv;
    json["invariant_point"]["CART"] = op.get_location(CART);
    json["invariant_point"]["FRAC"] = op.get_location(FRAC);

    // enum symmetry_type {identity_op, mirror_op, glide_op, rotation_op, screw_op, inversion_op, rotoinversion_op, invalid_op};
    if(op.type() == SymOp::identity_op) {
      json["type"] = "identity";
    }
    else if(op.type() == SymOp::mirror_op) {
      json["type"] = "mirror";
      json["mirror_normal"]["CART"] = op.get_eigenvec(CART);
      json["mirror_normal"]["FRAC"] = op.get_eigenvec(FRAC);
    }
    else if(op.type() == SymOp::glide_op) {
      json["type"] = "glide";
      json["mirror_normal"]["CART"] = op.get_eigenvec(CART);
      json["mirror_normal"]["FRAC"] = op.get_eigenvec(FRAC);
      json["shift"]["CART"] = op.get_screw_glide_shift(CART);
      json["shift"]["FRAC"] = op.get_screw_glide_shift(FRAC);
    }
    else if(op.type() == SymOp::rotation_op) {
      json["type"] = "rotation";
      json["rotation_axis"]["CART"] = op.get_eigenvec(CART);
      json["rotation_axis"]["FRAC"] = op.get_eigenvec(FRAC);
      json["rotation_angle"] = op.get_rotation_angle();
    }
    else if(op.type() == SymOp::screw_op) {
      json["type"] = "screw";
      json["rotation_axis"]["CART"] = op.get_eigenvec(CART);
      json["rotation_axis"]["FRAC"] = op.get_eigenvec(FRAC);
      json["rotation_angle"] = op.get_rotation_angle();
      json["shift"]["CART"] = op.get_screw_glide_shift(CART);
      json["shift"]["FRAC"] = op.get_screw_glide_shift(FRAC);
    }
    else if(op.type() == SymOp::rotoinversion_op) {
      json["type"] = "rotoinversion";
      json["rotation_axis"]["CART"] = op.get_eigenvec(CART);
      json["rotation_axis"]["FRAC"] = op.get_eigenvec(FRAC);
      json["rotation_angle"] = op.get_rotation_angle();
    }
    else if(op.type() == SymOp::invalid_op) {
      json["type"] = "invalid";
    }

  }

  inline void write_symgroup(const SymGroup &grp, jsonParser &json) {
    json = jsonParser::object();

    json["symop"] = jsonParser::array(grp.size());
    for(int i = 0; i < grp.size(); i++) {
      write_symop(grp[i], json["symop"][i], grp.class_of_op(i), grp.ind_inverse(i));
    }
    json["name"] = grp.get_name();
    json["latex_name"] = grp.get_latex_name();
    json["periodicity"] = grp.get_periodicity();
    if(grp.get_periodicity() == PERIODIC) {
      json["possible_space_groups"] = grp.possible_space_groups();
    }
    json["conjugacy_class"] = grp.get_conjugacy_classes();
    json["inverse"] = jsonParser::array(grp.size());
    for(int i = 0; i < grp.size(); i++) {
      json["inverse"][i] = grp.ind_inverse(i);
    }
    json["multiplication_table"] = grp.get_multi_table();
    json["character_table"] = grp.get_character_table();
  }

  
  // --------- ChemicalReference IO Declarations --------------------------------------------------
  
  ChemicalReference read_chemical_reference(fs::path filename, const Structure& prim, double tol = 1e-14);

  ChemicalReference read_chemical_reference(const jsonParser& json, const Structure& prim, double tol = 1e-14);
  
  void write_chemical_reference(const ChemicalReference& chem_ref, fs::path filename);

  void write_chemical_reference(const ChemicalReference& chem_ref, const jsonParser& json);
  
  
  // --------- ChemicalReference IO Definitions --------------------------------------------------
  
  /// \brief Read chemical reference states from JSON file
  ///
  /// See documentation in related function for expected form of the JSON
  inline ChemicalReference read_chemical_reference(fs::path filename,
                                            const Structure& prim, 
                                            double tol) {
    try {
      jsonParser json(filename);
      return read_chemical_reference(json, prim, tol);
    }
    catch(...) {
      std::cerr << "Error reading chemical reference states from " << filename << std::endl;
      /// re-throw exceptions
      throw;
    }
  }
  
  /// \brief Read chemical reference states from JSON
  ///
  /// Example expected form:
  /// \code
  /// { 
  ///   "chemical_reference" : {
  ///     "global" : ...,
  ///     "supercell": {
  ///       "SCELX": ...,
  ///       "SCELY": ...
  ///     },
  ///     "config": {
  ///       "SCELX/I": ...,
  ///       "SCELY/J": ...
  ///     }
  ///   }
  /// }
  /// \endcode
  ///
  /// See one_chemical_reference_from_json for documentation of the \code {...} \endcode expected form.
  inline ChemicalReference read_chemical_reference(const jsonParser& json,
                                            const Structure& prim, 
                                            double tol) {
    
    if(json.find("chemical_reference") == json.end()) {
      throw std::runtime_error("Error reading chemical reference states: Expected \"chemical_reference\" entry");
    }
    
    return jsonConstructor<ChemicalReference>::from_json(json["chemical_reference"], prim, tol);
    
  }
  
  inline void write_chemical_reference(const ChemicalReference& chem_ref, fs::path filename) {
    SafeOfstream outfile;
    outfile.open(filename);

    jsonParser json;
    write_chemical_reference(chem_ref, json);
    json.print(outfile.ofstream());

    outfile.close();
  }

  inline void write_chemical_reference(const ChemicalReference& chem_ref, jsonParser& json) {
    json.put_obj();
    to_json(chem_ref, json["chemical_reference"]);
  }


  // --------- CompositionAxes Declarations --------------------------------------------------

  struct CompositionAxes {

    CompositionAxes() {}

    /// \brief Read CompositionAxes from file
    CompositionAxes(fs::path filename);

    /// \brief Read CompositionAxes from JSON
    CompositionAxes(const jsonParser &json);

    /// \brief Read CompositionAxes from file
    void read(fs::path filename);

    /// \brief Read CompositionAxes from JSON
    void read(const jsonParser &json, fs::path filename = fs::path());

    /// \brief Write CompositionAxes to file
    void write(fs::path filename) const;

    /// \brief Write CompositionAxes to JSON
    void write(jsonParser &json) const;


    std::map<std::string, CompositionConverter> standard;
    std::map<std::string, CompositionConverter> custom;
    bool has_current_axes = false;
    std::string curr_key;
    CompositionConverter curr;

    int err_code = 0;
    std::string err_message;
  };


  /// \brief Read standard axes from JSON, and output to std::map<std::string, CompositionConverter>
  template<typename OutputIterator>
  OutputIterator read_composition_axes(OutputIterator result, const jsonParser &json);


  // --------- CompositionAxes Definitions --------------------------------------------------

  /// \brief Read CompositionAxes from file
  ///
  /// - Stores err_code and err_message for non-fatal errors
  /// - Catches and throws any exceptions
  inline CompositionAxes::CompositionAxes(fs::path filename) {

    try {

      read(jsonParser(filename), filename);

    }
    catch(...) {
      std::cerr << "Error reading composition axes from " << filename << std::endl;
      throw;
    }
  }

  /// \brief Read CompositionAxes from file
  ///
  /// - Stores err_code and err_message for non-fatal errors
  /// - Catches and throws any exceptions
  inline CompositionAxes::CompositionAxes(const jsonParser &json) {
    read(json);
  }

  inline void CompositionAxes::read(fs::path filename) {

    try {

      read(jsonParser(filename), filename);

    }
    catch(...) {

      std::cerr << "Error reading composition axes from " << filename << std::endl;
      throw;
    }
  }

  inline void CompositionAxes::read(const jsonParser &json, fs::path filename) {

    try {

      *this = CompositionAxes();

      if(json.contains("standard_axes")) {
        read_composition_axes(std::inserter(standard, standard.begin()), json["standard_axes"]);
      }

      if(json.contains("custom_axes")) {
        read_composition_axes(std::inserter(custom, custom.begin()), json["custom_axes"]);
      }

      has_current_axes = json.get_if(curr_key, "current_axes");

      if(has_current_axes) {

        if(standard.find(curr_key) != standard.cend() &&
           custom.find(curr_key) != custom.cend()) {

          std::string tmp = filename.empty() ? "the JSON" : filename.string();

          std::stringstream ss;
          ss << "Error: The current composition axes specified in " << tmp <<
             " can found in both the standard and custom compostion axes.\n\n" <<
             "Please edit the custom composition axes to remove this ambiguity.";

          err_message = ss.str();

          err_code = 1;

        }
        else if(standard.find(curr_key) == standard.cend() &&
                custom.find(curr_key) == custom.cend()) {

          std::string tmp = filename.empty() ? "the JSON" : filename.string();

          std::stringstream ss;
          ss << "Warning: The current composition axes specified in " << tmp <<
             " can not be found in the standard or custom compostion axes.\n\n" <<
             "Please use 'casm composition --select' to re-select your composition axes,\n" <<
             "Please use 'casm composition --calc' to re-calc your standard axes,\n" <<
             "or edit the custom composition axes.";

          err_message = ss.str();

          err_code = 2;
        }
        else if(standard.find(curr_key) != standard.cend()) {
          curr = standard[curr_key];
        }
        else if(custom.find(curr_key) != custom.cend()) {
          curr = custom[curr_key];
        }
      }

    }
    catch(...) {
      throw;
    }
  }

  /// \brief Write CompositionAxes to file
  inline void CompositionAxes::write(fs::path filename) const {

    SafeOfstream outfile;
    outfile.open(filename);
    jsonParser json;
    write(json);
    json.print(outfile.ofstream());
    outfile.close();
  }

  /// \brief Write CompositionAxes to JSON
  inline void CompositionAxes::write(jsonParser &json) const {

    json = jsonParser::object();

    json["standard_axes"][standard.cbegin()->first] = true;

    for(auto it = standard.cbegin(); it != standard.cend(); ++it) {
      json["standard_axes"][it->first] = it->second;
    }

    for(auto it = custom.cbegin(); it != custom.cend(); ++it) {
      json["custom_axes"][it->first] = it->second;
    }

    if(has_current_axes) {
      json["current_axes"] = curr_key;
    }

  }

  /// \brief Read standard axes from JSON, and output to std::map<std::string, CompositionConverter>
  ///
  /// - json Format: {"standard_axes": {"key0" : {CompositionConverter}, ... more axes ...},
  ///                 "custom_axes":   {"key0" : {CompositionConverter}, ... more axes ...},
  ///                 "current_axes": "key"}
  ///   - "axes_type" is either "standard_axes" or "custom_axes"
  ///   - "key" is a string indicating which composition axes in the "standard_axes" or "custom_axes" JSON object
  ///
  /// - This read json["standard_axes"] or json["custom_axes"]
  ///
  template<typename OutputIterator>
  OutputIterator read_composition_axes(OutputIterator result, const jsonParser &json) {

    CompositionConverter conv;
    for(auto it = json.cbegin(); it != json.cend(); ++it) {
      from_json(conv, *it);
      *result++ = std::make_pair(it.name(), conv);
    }
    return result;
  }


  // ---------- prim_nlist.json IO -------------------------------------------------------------

  inline void write_prim_nlist(const Array<UnitCellCoord> &prim_nlist, const fs::path &nlistpath) {

    SafeOfstream outfile;
    outfile.open(nlistpath);
    jsonParser nlist_json;
    nlist_json = prim_nlist;
    nlist_json.print(outfile.ofstream());
    outfile.close();

    return;
  }

  inline Array<UnitCellCoord> read_prim_nlist(const fs::path &nlistpath) {

    try {
      Array<UnitCellCoord> prim_nlist;
      jsonParser read_json(nlistpath);
      from_json(prim_nlist, read_json);
      return prim_nlist;
    }
    catch(...) {
      std::cerr << "Error reading: " << nlistpath;
      throw;
    }
  }



}

#endif
