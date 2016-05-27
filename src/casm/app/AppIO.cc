#include "casm/app/AppIO.hh"
#include "casm/symmetry/SymInfo.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/clusterography/jsonClust.hh"

namespace CASM {


  // --------- PrimIO Definitions --------------------------------------------------

  BasicStructure<Site> read_prim(fs::path filename) {

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
  BasicStructure<Site> read_prim(const jsonParser &json) {

    try {

      // read lattice
      Eigen::Matrix3d latvec_transpose;

      from_json(latvec_transpose, json["lattice_vectors"]);

      Lattice lat(latvec_transpose.transpose());

      // create prim using lat
      BasicStructure<Site> prim(lat);

      // read title
      from_json(prim.title, json["title"]);

      Eigen::Vector3d vec;

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
      for(int i = 0; i < json["basis"].size(); i++) {

        // read coordinate
        Eigen::Vector3d coord(json["basis"][i]["coordinate"][0].get<double>(),
                              json["basis"][i]["coordinate"][1].get<double>(),
                              json["basis"][i]["coordinate"][2].get<double>());
        Site site(prim.lattice());
        if(mode == FRAC) {
          site.frac() = coord;
        }
        else if(mode == CART) {
          site.cart() = coord;
        }

        // read atom occupant names
        Array<std::string> occ_name;
        from_json(occ_name, json["basis"][i]["occupant_dof"]);

        // fill site.site_occupant
        Array<Molecule> tocc;
        for(int i = 0; i < occ_name.size(); i++) {
          Molecule tMol(prim.lattice());
          tMol.name = occ_name[i];
          tMol.push_back(AtomPosition(0, 0, 0, occ_name[i], prim.lattice(), CART));
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
  void write_prim(const BasicStructure<Site> &prim, fs::path filename, COORD_TYPE mode) {

    SafeOfstream outfile;
    outfile.open(filename);

    jsonParser json;
    write_prim(prim, json, mode);
    json.print(outfile.ofstream());

    outfile.close();
  }

  /// \brief Write prim.json as JSON
  void write_prim(const BasicStructure<Site> &prim, jsonParser &json, COORD_TYPE mode) {

    json = jsonParser::object();

    json["title"] = prim.title;

    json["lattice_vectors"] = prim.lattice().lat_column_mat().transpose();

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
      json["basis"][i]["coordinate"].put_array();
      if(mode == FRAC) {
        json["basis"][i]["coordinate"].push_back(prim.basis[i].frac(0));
        json["basis"][i]["coordinate"].push_back(prim.basis[i].frac(1));
        json["basis"][i]["coordinate"].push_back(prim.basis[i].frac(2));
      }
      else if(mode == CART) {
        json["basis"][i]["coordinate"].push_back(prim.basis[i].cart(0));
        json["basis"][i]["coordinate"].push_back(prim.basis[i].cart(1));
        json["basis"][i]["coordinate"].push_back(prim.basis[i].cart(2));
      }

      json["basis"][i]["occupant_dof"] = jsonParser::array(prim.basis[i].site_occupant().size());

      for(int j = 0; j < prim.basis[i].site_occupant().size(); j++) {
        json["basis"][i]["occupant_dof"][j] = prim.basis[i].site_occupant()[j].name;
      }

    }
  }


  // --------- SymmetryIO Declarations --------------------------------------------------

  void write_symop(const SymGroup &grp, Index i, jsonParser &j) {
    j = jsonParser::object();

    const SymOp &op = grp[i];

    to_json(op.matrix(), j["matrix"]["CART"]);
    to_json(grp.lattice().inv_lat_column_mat()*op.matrix()*grp.lattice().lat_column_mat(), j["matrix"]["FRAC"]);

    to_json_array(op.tau(), j["tau"]["CART"]);
    to_json_array(grp.lattice().inv_lat_column_mat()*op.tau(), j["tau"]["FRAC"]);

    to_json(grp.class_of_op(i), j["conjugacy_class"]);
    to_json(grp.ind_inverse(i), j["inverse"]);

    add_sym_info(grp.info(i), j);
  }

  void write_symgroup(const SymGroup &grp, jsonParser &json) {
    json = jsonParser::object();

    json["symop"] = jsonParser::array(grp.size());
    for(int i = 0; i < grp.size(); i++) {
      write_symop(grp, i, json["symop"][i]);
    }

    json["name"] = grp.get_name();
    json["latex_name"] = grp.get_latex_name();
    json["periodicity"] = grp.periodicity();
    if(grp.periodicity() == PERIODIC) {
      json["possible_space_groups"] = grp.possible_space_groups();
    }
    json["conjugacy_class"] = grp.get_conjugacy_classes();
    json["inverse"] = jsonParser::array(grp.size());
    for(int i = 0; i < grp.size(); i++) {
      json["inverse"][i] = grp.ind_inverse(i);
    }
    json["multiplication_table"] = grp.get_multi_table();
    json["character_table"] = grp.character_table();
  }


  // --------- ChemicalReference IO Definitions --------------------------------------------------

  /// \brief Read chemical reference states from JSON file
  ///
  /// See documentation in related function for expected form of the JSON
  ChemicalReference read_chemical_reference(fs::path filename,
                                            const Structure &prim,
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
  ChemicalReference read_chemical_reference(const jsonParser &json,
                                            const Structure &prim,
                                            double tol) {

    if(json.find("chemical_reference") == json.end()) {
      throw std::runtime_error("Error reading chemical reference states: Expected \"chemical_reference\" entry");
    }

    return jsonConstructor<ChemicalReference>::from_json(json["chemical_reference"], prim, tol);

  }

  void write_chemical_reference(const ChemicalReference &chem_ref, fs::path filename) {
    SafeOfstream outfile;
    outfile.open(filename);

    jsonParser json;
    write_chemical_reference(chem_ref, json);
    json.print(outfile.ofstream());

    outfile.close();
  }

  void write_chemical_reference(const ChemicalReference &chem_ref, jsonParser &json) {
    json.put_obj();
    to_json(chem_ref, json["chemical_reference"]);
  }


  // --------- CompositionAxes Definitions --------------------------------------------------

  /// \brief Read CompositionAxes from file
  ///
  /// - Stores err_code and err_message for non-fatal errors
  /// - Catches and throws any exceptions
  CompositionAxes::CompositionAxes(fs::path _filename) {
    read(_filename);
  }

  /// \brief Read CompositionAxes from file
  ///
  /// - Stores err_code and err_message for non-fatal errors
  /// - Catches and throws any exceptions
  CompositionAxes::CompositionAxes(const jsonParser &json) {
    read(json);
  }

  void CompositionAxes::read(fs::path _filename) {

    try {

      read(jsonParser(_filename), _filename);

    }
    catch(...) {

      std::cerr << "Error reading composition axes from " << _filename << std::endl;
      throw;
    }
  }

  void CompositionAxes::read(const jsonParser &json, fs::path _filename) {

    try {

      *this = CompositionAxes();

      filename = _filename;

      filename = _filename;

      if(json.contains("standard_axes")) {
        read_composition_axes(std::inserter(standard, standard.begin()), json["standard_axes"]);
      }

      if(json.contains("custom_axes")) {
        read_composition_axes(std::inserter(custom, custom.begin()), json["custom_axes"]);
      }

      std::string key;
      has_current_axes = json.get_if(key, "current_axes");

      if(has_current_axes) {
        select(key);
      }

    }
    catch(...) {
      throw;
    }
  }

  /// \brief Write CompositionAxes to file
  void CompositionAxes::write(fs::path _filename) {

    filename = _filename;
    SafeOfstream outfile;
    outfile.open(filename);
    jsonParser json;
    write(json);
    json.print(outfile.ofstream());
    outfile.close();
  }

  /// \brief Write CompositionAxes to JSON
  void CompositionAxes::write(jsonParser &json) const {

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

  /// \brief Set this->curr using key
  void CompositionAxes::select(std::string key) {
    if(standard.find(key) != standard.cend() &&
       custom.find(key) != custom.cend()) {

      std::stringstream ss;
      ss << "Error: The composition axes " << key <<
         " can found in both the standard and custom compostion axes.\n\n" <<
         "Please edit the custom composition axes to remove this ambiguity.";

      err_message = ss.str();

      err_code = 1;

    }
    else if(standard.find(key) == standard.cend() &&
            custom.find(key) == custom.cend()) {

      std::stringstream ss;
      ss << "Warning: The composition axes " << key <<
         " can not be found in the standard or custom compostion axes.\n\n" <<
         "Please use 'casm composition --select' to re-select your composition axes,\n" <<
         "Please use 'casm composition --calc' to re-calc your standard axes,\n" <<
         "or edit the custom composition axes.";

      err_message = ss.str();

      err_code = 2;
    }
    else if(standard.find(key) != standard.cend()) {
      curr = standard[key];
    }
    else if(custom.find(key) != custom.cend()) {
      curr = custom[key];
    }
    curr_key = key;
    has_current_axes = true;

  }

  // ---------- prim_nlist.json IO -------------------------------------------------------------

  void write_prim_nlist(const Array<UnitCellCoord> &prim_nlist, const fs::path &nlistpath) {

    SafeOfstream outfile;
    outfile.open(nlistpath);
    jsonParser nlist_json;
    nlist_json = prim_nlist;
    nlist_json.print(outfile.ofstream());
    outfile.close();

    return;
  }

  Array<UnitCellCoord> read_prim_nlist(const fs::path &nlistpath) {

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


  // ---------- basis.json IO ------------------------------------------------------------------

  /// \brief Write summary of basis functions
  ///
  /// Format:
  /// \code
  /// {
  ///   "site_functions":[
  ///     {
  ///       "asym_unit": X,
  ///       "sublat_indices: [2, 3],
  ///       "phi_b_0": {"Va":0.0, "O":1.0},
  ///       "phi_b_1": {"Va":0.0, "O":1.0},
  ///        ...
  ///     },
  ///     ...
  ///   ],
  ///   "cluster_functions":[
  ///     {
  ///       "eci": X.XXXXX,
  ///       "prototype_function": "\phi_b_i(s_j)...",
  ///       "orbit": [branch_index, orbit_index],
  ///       "linear_orbit_index": I,
  ///       "mult": X,
  ///       "prototype": [
  ///         [b, i, j, k],
  ///         ...
  ///       ]
  ///     },
  ///     ...
  ///   ]
  /// }
  /// \endcode
  ///
  void write_basis(const SiteOrbitree &tree, const Structure &prim, jsonParser &json, double tol) {

    json = jsonParser::object();

    //   "site_functions":[
    //     {
    //       "asym_unit": X,
    //       "sublat": [2, 3],
    //       "basis": {
    //         "phi_b_0": {"Va":0.0, "O":1.0},
    //         "phi_b_1": {"Va":0.0, "O":1.0}
    //       }
    //     },
    //     ...
    //   ],

    jsonParser &sitef = json["site_functions"];
    sitef = jsonParser::array(prim.basis.size(), jsonParser::object());
    for(Index no = 0; no < tree.asym_unit().size(); no++) {
      for(Index ne = 0; ne < tree.asym_unit()[no].size(); ne++) {

        const SiteCluster &equiv = tree.asym_unit()[no][ne];
        const Site &site = equiv[0];

        Index b = site.basis_ind();
        sitef[b]["sublat"] = b;
        sitef[b]["asym_unit"] = no;

        if(equiv.clust_basis.size() == 0) {
          sitef[b]["basis"].put_null();
        }
        else {
          for(Index f = 0; f < equiv.clust_basis.size(); f++) {
            std::stringstream fname;
            fname << "\\phi_" << b << '_' << f;
            for(Index s = 0; s < site.site_occupant().size(); s++) {

              // "\phi_b_f": {"Zr":0.0, ...}
              sitef[b]["basis"][fname.str()][site.site_occupant()[s].name] =
                equiv.clust_basis[f]->eval(
                  Array<Index>(1, site.site_occupant().ID()),
                  Array<Index>(1, s)
                );
            }
          }
        }
      }
    }

    //   "cluster_functions":[
    //     {
    //       ("eci": X.XXXXX,) <-- is included after fitting
    //       "prototype_function": "\phi_b_i(s_j)...",
    //       "orbit": [branch_index, cluster_orbit_index, bfunc_index],
    //       "linear_orbit_index": I,
    //       "mult": X,
    //       "prototype": {
    //         "max_length": X.X,
    //         "min_length": X.X,
    //         "sites": [
    //           [b, i, j, k],
    //           ...
    //         ]
    //       }
    //     },
    //     ...
    //   ]
    // }

    jsonParser &orbitf = json["cluster_functions"];
    orbitf = jsonParser::array();
    for(Index i = 0; i < tree.size(); i++) {
      jsonParser tjson;
      for(Index j = 0; j < tree.size(i); j++) { //Loops over all i sized Orbits of clusters

        // cluster specific info
        tjson["orbit"] = std::vector<Index>({i, j, 0});
        tjson["mult"] = tree.orbit(i, j).size();
        to_json(jsonHelper(tree.orbit(i, j)[0], prim), tjson["prototype"]);

        // basis function info
        BasisSet tbasis(tree.orbit(i, j)[0].clust_basis);
        tbasis.accept(OccFuncLabeler("\\phi_%b_%f(s_%n)"));

        for(Index nf = 0; nf < tbasis.size(); ++nf) {
          tjson["orbit"][2] = nf;
          tjson["prototype_function"] = tbasis[nf]->tex_formula();
          orbitf.push_back(tjson);
        }
      }
    }

    for(Index i = 0; i < orbitf.size(); ++i) {
      orbitf[i]["linear_function_index"] = i;
    }
  }

}

