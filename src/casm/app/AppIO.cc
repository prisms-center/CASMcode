#include "casm/app/AppIO_impl.hh"
#include "casm/symmetry/SymInfo.hh"
#include "casm/basis_set/FunctionVisitor.hh"

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
  /*
    void write_prim_nlist(const Array<UnitCellCoord> &prim_nlist, const fs::path &nlistpath) {

      SafeOfstream outfile;
      outfile.open(nlistpath);
      jsonParser nlist_json;
      nlist_json = prim_nlist;
      nlist_json.print(outfile.ofstream());
      outfile.close();

      return;
    }

    Array<UnitCellCoord> read_prim_nlist(const fs::path &nlistpath, const BasicStructure<Site>& prim) {

      try {
        jsonParser json(nlistpath);
        Array<UnitCellCoord> prim_nlist(json.size(), UnitCellCoord(prim));
        auto nlist_it = prim_nlist.begin();
        for(auto it=json.begin(); it!=json.end(); ++it) {
          from_json(*nlist_it++, *it);
        }
        return prim_nlist;
      }
      catch(...) {
        std::cerr << "Error reading: " << nlistpath;
        throw;
      }
    }
  */

  // ---------- Orbit<IntegralCluster> & ClexBasis IO ------------------------------------------------------------------


  SitesPrinter::SitesPrinter(int _indent_space, char _delim, COORD_TYPE _mode) :
    indent_space(_indent_space),
    delim(_delim),
    mode(_mode) {}

  std::string SitesPrinter::indent() const {
    return std::string(indent_space, ' ');
  }

  void SitesPrinter::coord_mode(std::ostream &out) {
    out << "COORD_MODE = " << mode << std::endl << std::endl;
  }

  void SitesPrinter::print_sites(const IntegralCluster &clust, std::ostream &out) {
    for(const auto &coord : clust) {
      out << indent() << indent() << indent();
      out.setf(std::ios::showpoint, std::ios_base::fixed);
      out.precision(5);
      out.width(9);
      coord.site().print(out);
      if(delim)
        out << delim;
      out << std::flush;
    }
  }


  // explicit template instantiations

#define PRINT_CLUST_INST(ITERATOR,INSERTER,PRINTER) \
  template void print_clust<ITERATOR,PRINTER>(ITERATOR begin, ITERATOR end, std::ostream &out, PRINTER printer); \
  template jsonParser &write_clust<ITERATOR>(ITERATOR begin, ITERATOR end, jsonParser &json, PRINTER printer); \
  template jsonParser &write_clust<ITERATOR>(ITERATOR begin, ITERATOR end, jsonParser &json, PRINTER printer, const jsonParser &bspecs);

#define ORBIT_CONTAINER_INST(ITERATOR,INSERTER,ORBIT) \
  PRINT_CLUST_INST(ITERATOR,INSERTER,ProtoSitesPrinter) \
  PRINT_CLUST_INST(ITERATOR,INSERTER,FullSitesPrinter) \
  PRINT_CLUST_INST(ITERATOR,INSERTER,ProtoFuncsPrinter) \
  template void print_site_basis_funcs<ITERATOR>(ITERATOR begin, ITERATOR end, const ClexBasis &clex_basis, std::ostream &out, COORD_TYPE mode); \
  template INSERTER read_clust<INSERTER, typename ORBIT::SymCompareType>(INSERTER result, const jsonParser &json, const Structure &prim, const SymGroup& generating_grp, const typename ORBIT::SymCompareType &sym_compare);

#define _VECTOR_IT(ORBIT) std::vector<ORBIT>::iterator
#define _VECTOR_INSERTER(ORBIT) std::back_insert_iterator<std::vector<ORBIT> >

#define _SET_IT(ORBIT) std::set<ORBIT>::iterator
#define _SET_INSERTER(ORBIT) std::insert_iterator<std::set<ORBIT> >

#define ORBIT_VECTOR_INST(ORBIT) ORBIT_CONTAINER_INST(_VECTOR_IT(ORBIT),_VECTOR_INSERTER(ORBIT), ORBIT)
#define ORBIT_SET_INST(ORBIT) ORBIT_CONTAINER_INST(_SET_IT(ORBIT),_SET_INSERTER(ORBIT), ORBIT)

  ORBIT_VECTOR_INST(LocalIntegralClusterOrbit)
  ORBIT_VECTOR_INST(PrimPeriodicIntegralClusterOrbit)
  ORBIT_VECTOR_INST(ScelPeriodicIntegralClusterOrbit)

  ORBIT_SET_INST(LocalIntegralClusterOrbit)
  ORBIT_SET_INST(PrimPeriodicIntegralClusterOrbit)
  ORBIT_SET_INST(ScelPeriodicIntegralClusterOrbit)

}

