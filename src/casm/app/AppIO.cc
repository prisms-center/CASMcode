#include "casm/app/AppIO_impl.hh"

#include "casm/app/HamiltonianModules.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/crystallography/BasicStructure_impl.hh"
#include "casm/clex/io/json/ChemicalReference.hh"
#include "casm/clusterography/ClusterSymCompare.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/global/enum/json_io.hh"
#include "casm/global/enum/stream_io.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/symmetry/SymInfo.hh"
#include "casm/symmetry/Orbit_impl.hh"

namespace CASM {

  using xtal::SpeciesAttribute;

  // --------- PrimIO Definitions --------------------------------------------------

  jsonParser const &from_json(SpeciesAttribute &_attr, jsonParser const &json) {
    _attr.set_value(json["value"].get<Eigen::VectorXd>());
    return json;
  }

  //****************************************************

  jsonParser &to_json(SpeciesAttribute const &_attr, jsonParser &json) {
    json.put_obj();
    to_json_array(_attr.value(), json["value"]);
    return json;
  }

  //****************************************************
  jsonParser &to_json(AtomPosition const &apos, jsonParser &json, Eigen::Ref<const Eigen::Matrix3d> const &c2f_mat) {
    json.put_obj();
    to_json_array(c2f_mat * apos.cart(), json["coordinate"]);
    json["name"] = apos.name();
    if(apos.attributes().size())
      json["attributes"] = apos.attributes();
    return json;
  }

  //****************************************************

  void from_json(AtomPosition &apos, const jsonParser &json,  Eigen::Ref<const Eigen::Matrix3d> const &f2c_mat, HamiltonianModules const &_modules) {
    apos = json.get<AtomPosition>(f2c_mat, _modules);

    return;
  }

  //****************************************************

  AtomPosition jsonConstructor<AtomPosition>::from_json(const jsonParser &json, Eigen::Matrix3d const &f2c_mat, HamiltonianModules const &_modules) {
    std::string _name;
    Eigen::Vector3d _pos(0., 0., 0.);
    std::map<std::string, SpeciesAttribute> attr_map;
    if(json.is_obj()) {
      _name = json["name"].get<std::string>();
      if(json.contains("coordinate")) {
        _pos = f2c_mat * json["coordinate"].get<Eigen::Vector3d>();
        //std::cout << "f2c_mat: \n" << f2c_mat << "\n";
        //std::cout << "_pos: " << _pos.transpose() << "\n";
      }
      if(json.contains("attributes")) {
        auto it = json["attributes"].cbegin(), end_it = json["attributes"].cend();
        for(; it != end_it; ++it) {
          auto result_pair = attr_map.emplace(it.name(), _modules.aniso_val_dict().lookup(it.name()));
          CASM::from_json(result_pair.first->second, *it);
        }
      }

    }
    else if(json.is_string()) {
      _name = json.get<std::string>();
    }
    else
      throw std::runtime_error("Invalid JSON input encountered. Unable to parse AtomPosition object.\n");

    AtomPosition result(_pos, _name);
    result.set_attributes(attr_map);
    return result;
  }

  //****************************************************
  //   Write Molecule to json.
  //****************************************************

  jsonParser &to_json(Molecule const &mol, jsonParser &json, Eigen::Ref<const Eigen::Matrix3d> const &c2f_mat) {
    json.put_obj();
    CASM::to_json(mol.atoms(), json["atoms"], c2f_mat);
    json["name"] = mol.name();
    if(mol.attributes().size())
      json["attributes"] = mol.attributes();

    return json;
  }

  //****************************************************
  //
  //    Read Molecule from json.
  //
  //****************************************************

  void from_json(Molecule &mol, const jsonParser &json, Eigen::Ref<const Eigen::Matrix3d> const &f2c_mat, HamiltonianModules const &_modules) {
    std::vector<AtomPosition> _atoms;
    if(json.contains("atoms")) {
      CASM::from_json(_atoms, json["atoms"], f2c_mat, _modules);
    }
    mol.set_atoms(_atoms);

    std::map<std::string, SpeciesAttribute> attr_map;
    if(json.contains("attributes")) {
      auto it = json["attributes"].cbegin(), end_it = json["attributes"].cend();
      for(; it != end_it; ++it) {
        auto result_pair = attr_map.emplace(it.name(), _modules.aniso_val_dict().lookup(it.name()));
        from_json(result_pair.first->second, *it);
      }
    }
    mol.set_attributes(attr_map);

    //jsonParser tjson;
    //to_json(mol,tjson,f2c_mat.inverse());
    //std::cout << "Read Molecule :\n"<< json << "\n";
  }

  //****************************************************
  //
  //****************************************************

  Molecule jsonConstructor<Molecule>::from_json(const jsonParser &json,  Eigen::Ref<const Eigen::Matrix3d> const &f2c_mat, HamiltonianModules const &_modules) {
    return json.get<Molecule>(f2c_mat, _modules);
  }

  //****************************************************

  Site jsonConstructor<Site>::from_json(const jsonParser &json,
                                        Lattice const &_home,
                                        COORD_TYPE coordtype,
                                        std::map<std::string, Molecule> const &mol_map,
                                        HamiltonianModules const &_modules) {
    Site result(_home);
    CASM::from_json(result, json, _home, coordtype, mol_map, _modules);
    return result;
  }

  //****************************************************

  jsonParser &to_json(const Site &site, jsonParser &json, COORD_TYPE coordtype) {
    json.put_obj();

    // class Site : public Coordinate
    if(coordtype == FRAC)
      to_json_array(site.frac(), json["coordinate"]);
    else
      to_json_array(site.cart(), json["coordinate"]);

    // Occupant_DoF<Molecule> occupant_dof;
    Eigen::Matrix3d c2f = Eigen::Matrix3d::Identity();
    // change this to use FormatFlag
    if(coordtype == FRAC)
      c2f = site.home().inv_lat_column_mat();
    CASM::to_json(site.occupant_dof().domain(), json["occupants"], c2f);

    if(site.dofs().size()) {
      json["dofs"] = site.dofs();
    }

    // Index m_label
    if(valid_index(site.label()))
      json["label"] = site.label();


    return json;

  }

  //****************************************************

  void from_json(Site &site,
                 const jsonParser &json,
                 Lattice const &_home,
                 COORD_TYPE coordtype,
                 std::map<std::string, Molecule> const &mol_map,
                 HamiltonianModules const &_modules) {
    site.set_lattice(_home, coordtype);
    if(coordtype == FRAC)
      site.frac() = json["coordinate"].get<Eigen::Vector3d>();
    else
      site.cart() = json["coordinate"].get<Eigen::Vector3d>();

    // Index m_label -- must be greater than zero
    Index _label = -1;
    if(json.contains("label")) {
      CASM::from_json(_label, json["label"]);
      if(!valid_index(_label))
        throw std::runtime_error("JSON specification of site has {\"label\" : " + std::to_string(_label) + "}, but \"label\" must be greater than 0.\n");
    }
    site.set_label(_label);

    // Local continuous dofs

    std::map<std::string, DoFSet> _dof_map;
    if(json.contains("dofs")) {

      auto it = json["dofs"].begin(), end_it = json["dofs"].end();
      for(; it != end_it; ++it) {
        if(_dof_map.count(it.name()))
          throw std::runtime_error("Error parsing global field \"dofs\" from JSON. DoF type " + it.name() + " cannot be repeated.");

        try {
          _dof_map.emplace(std::make_pair(it.name(), it->get<DoFSet>(_modules.aniso_val_dict().lookup(it.name()))));
        }
        catch(std::exception &e) {
          throw std::runtime_error("Error parsing global field \"dofs\" from JSON. Failure for DoF type " + it.name() + ": " + e.what());
        }

      }
    }
    site.set_dofs(_dof_map);

    std::vector<Molecule> t_occ;
    std::string occ_key;
    if(json.contains("occupants"))
      occ_key = "occupants";
    else if(json.contains("occupant_dof"))
      occ_key = "occupant_dof";

    if(!occ_key.empty()) {
      for(std::string const &occ : json[occ_key].get<std::vector<std::string> >()) {
        //std::cout << "CREATING OCCUPANT " << occ << "\n";
        // Have convenience options for attributes like magnetic moment, etc?
        auto it = mol_map.find(occ);
        if(it != mol_map.end())
          t_occ.push_back(it->second);
        else
          t_occ.push_back(Molecule::make_atom(occ));
      }
    }
    //std::cout << "t_occ.size() = " << t_occ.size() << "\n";
    if(t_occ.empty())
      t_occ.push_back(Molecule::make_unknown());
    site.set_allowed_occupants(t_occ);
  }

  BasicStructure<Site> read_prim(fs::path filename, double xtal_tol, HamiltonianModules const *_modules) {

    try {
      jsonParser json(filename);
      return read_prim(json, xtal_tol, _modules);
    }
    catch(...) {
      std::cerr << "Error reading prim from " << filename << std::endl;
      /// re-throw exceptions
      throw;
    }
  }

  /// \brief Read prim.json
  BasicStructure<Site> read_prim(const jsonParser &json, double xtal_tol, HamiltonianModules const *_modules) {
    HamiltonianModules default_module;
    if(_modules == nullptr)
      _modules = &default_module;

    try {

      // read lattice
      Eigen::Matrix3d latvec_transpose;

      from_json(latvec_transpose, json["lattice_vectors"]);

      Lattice lat(latvec_transpose.transpose(), xtal_tol);

      // create prim using lat
      BasicStructure<Site> prim(lat);

      // read title
      prim.set_title(json["title"].get<std::string>());

      Eigen::Vector3d vec;

      // Global DoFs
      {
        std::map<std::string, DoFSet> _dof_map;
        if(json.contains("dofs")) {
          auto it = json["dofs"].begin(), end_it = json["dofs"].end();
          for(; it != end_it; ++it) {
            if(_dof_map.count(it.name()))
              throw std::runtime_error("Error parsing global field \"dofs\" from JSON. DoF type " + it.name() + " cannot be repeated.");

            try {
              _dof_map.emplace(std::make_pair(it.name(), it->get<DoFSet>(_modules->aniso_val_dict().lookup(it.name()))));
            }
            catch(std::exception &e) {
              throw std::runtime_error("Error parsing global field \"dofs\" from JSON. Failure for DoF type " + it.name() + ": " + e.what());
            }

          }
          prim.set_global_dofs(_dof_map);
        }
      }

      // read basis coordinate mode
      COORD_TYPE mode;
      from_json(mode, json["coordinate_mode"]);

      // Molecules
      std::map<std::string, Molecule> mol_map;
      Eigen::Matrix3d f2c;
      if(mode == FRAC)
        f2c = lat.lat_column_mat();
      else
        f2c.setIdentity();

      if(json.contains("species")) {
        auto it = json["species"].begin();
        auto it_end = json["species"].end();
        for(; it != it_end; ++it) {
          std::string chem_name = it.name();
          it->get_if(chem_name, "name");
          //std::cout << "chem_name: " << chem_name << "\n";
          auto mol_it = mol_map.emplace(it.name(), Molecule(chem_name)).first;
          from_json(mol_it->second, *it, f2c, *_modules);
        }
      }


      // read basis sites
      for(jsonParser const &bjson : json["basis"])
        prim.push_back(bjson.get<Site>(prim.lattice(), mode, mol_map, *_modules));

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

    json["title"] = prim.title();

    json["lattice_vectors"] = prim.lattice().lat_column_mat().transpose();

    if(mode == COORD_DEFAULT) {
      mode = COORD_MODE::CHECK();
    }

    Eigen::Matrix3d c2f_mat;
    if(mode == FRAC) {
      c2f_mat = prim.lattice().inv_lat_column_mat();
      json["coordinate_mode"] = "Fractional";
    }
    else if(mode == CART) {
      c2f_mat.setIdentity();
      json["coordinate_mode"] = "Cartesian";
    }

    // Global DoFs
    for(auto const &_dof : prim.global_dofs()) {
      json["dofs"][_dof.first] = _dof.second;
    }
    auto mol_names = allowed_molecule_unique_names(prim);
    jsonParser &bjson = (json["basis"] = jsonParser::array(prim.basis().size()));
    for(int i = 0; i < prim.basis().size(); i++) {
      bjson[i] = jsonParser::object();
      jsonParser &cjson = bjson[i]["coordinate"].put_array();
      if(mode == FRAC) {
        cjson.push_back(prim.basis()[i].frac(0));
        cjson.push_back(prim.basis()[i].frac(1));
        cjson.push_back(prim.basis()[i].frac(2));
      }
      else if(mode == CART) {
        cjson.push_back(prim.basis()[i].cart(0));
        cjson.push_back(prim.basis()[i].cart(1));
        cjson.push_back(prim.basis()[i].cart(2));
      }

      if(prim.basis()[i].dofs().size()) {
        bjson[i]["dofs"] = prim.basis()[i].dofs();
      }


      jsonParser &ojson = (bjson[i]["occupants"] = jsonParser::array(prim.basis()[i].occupant_dof().size()));

      for(int j = 0; j < prim.basis()[i].occupant_dof().size(); j++) {
        ojson[j] = mol_names[i][j];
        if(prim.basis()[i].occupant_dof()[j].name() != mol_names[i][j]
           || !prim.basis()[i].occupant_dof()[j].is_atomic()) {
          jsonParser &ojson = to_json(prim.basis()[i].occupant_dof()[j], json["species"][mol_names[i][j]], c2f_mat);
          //if(prim.basis()[i].occupant_dof()[j].name()!=mol_names[i][j])
          //ojson["name"]=prim.basis()[i].occupant_dof()[j].name();
        }
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
    bool has_time_reversal = false;
    for(SymOp const &op : grp) {
      if(op.time_reversal()) {
        has_time_reversal = true;
        break;
      }
    }
    if(!has_time_reversal) {
      // For now, only print character table for nonmagnetic groups (character table routine does not work for magnetic groups)
      json["character_table"] = grp.character_table();
    }
  }


  // --------- ChemicalReference IO Definitions --------------------------------------------------

  /// \brief Read chemical reference states from JSON file
  ///
  /// See documentation in related function for expected form of the JSON
  ChemicalReference read_chemical_reference(fs::path filename,
                                            BasicStructure<Site> const &prim,
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
  ChemicalReference read_chemical_reference(jsonParser const &json,
                                            BasicStructure<Site> const &prim,
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

  /// \brief Serialize CompositionConverter to JSON
  jsonParser &to_json(const CompositionConverter &f, jsonParser &json) {

    json = jsonParser::object();
    json["components"] = f.components();
    json["independent_compositions"] = f.independent_compositions();
    json["origin"] = f.origin();
    for(int i = 0; i < f.independent_compositions(); i++) {
      json[CompositionConverter::comp_var(i)] = f.end_member(i);
    }
    json["mol_formula"] = f.mol_formula();
    json["param_formula"] = f.param_formula();

    return json;
  }

  /// \brief Deserialize CompositionConverter from JSON
  void from_json(CompositionConverter &f, const jsonParser &json) {

    try {

      std::vector<std::string> components;
      Eigen::VectorXd origin;

      int independent_compositions;

      from_json(components, json["components"]);
      from_json(origin, json["origin"]);

      from_json(independent_compositions, json["independent_compositions"]);
      Eigen::MatrixXd end_members(components.size(), independent_compositions);
      Eigen::VectorXd tvec;
      for(int i = 0; i < independent_compositions; i++) {
        from_json(tvec, json[CompositionConverter::comp_var(i)]);
        end_members.col(i) = tvec;
      }

      f = CompositionConverter(components.begin(), components.end(), origin, end_members);
      return;
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

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

      read(jsonParser(_filename));

    }
    catch(...) {

      std::cerr << "Error reading composition axes from " << _filename << std::endl;
      throw;
    }
  }

  void CompositionAxes::read(const jsonParser &json) {

    try {

      *this = CompositionAxes();

      json.get_if(all_axes, "possible_axes");

      json.get_if(enumerated, "enumerated");

      //Backward compatibility
      std::map<std::string, CompositionConverter> custom;
      json.get_if(custom, "custom_axes");
      all_axes.insert(custom.begin(), custom.end());

      if(json.contains("standard_axes")) {
        std::map<std::string, CompositionConverter> standard;
        from_json(standard, json["standard_axes"]);
        for(auto const &el : standard)
          enumerated.insert(el.first);
        standard.insert(standard.begin(), standard.end());
      }

      std::string key;
      json.get_if(key, "current_axes");
      if(!key.empty()) {
        select(key);
      }

    }
    catch(...) {
      throw;
    }
  }

  void CompositionAxes::erase_enumerated() {
    for(std::string const &el : enumerated) {
      all_axes.erase(el);
      if(curr_key == el)
        curr_key.clear();
    }
    enumerated.clear();
  }

  /// \brief Write CompositionAxes to file
  void CompositionAxes::write(fs::path _filename) const {

    SafeOfstream outfile;
    outfile.open(_filename);
    jsonParser json;
    write(json);
    json.print(outfile.ofstream());
    outfile.close();
  }

  /// \brief Write CompositionAxes to JSON
  void CompositionAxes::write(jsonParser &json) const {

    json = jsonParser::object();

    if(has_current_axes()) {
      json["current_axes"] = curr_key;
    }

    json["possible_axes"] = all_axes;

    if(enumerated.size())
      json["enumerated"] = enumerated;
  }

  /// \brief Set this->curr using key
  void CompositionAxes::select(std::string key) {
    auto it = all_axes.find(key);
    if(it == all_axes.end()) {
      std::stringstream ss;
      ss << "Warning: The composition axes " << key <<
         " cannot be found among the posible composition axes.\n\n" <<
         "Please use 'casm composition --select' to re-select your composition axes,\n" <<
         "use 'casm composition --calc' to re-calc your standard axes,\n" <<
         "or add custom composition axes manually.";

      err_message = ss.str();

      err_code = 2;
    }
    else {
      curr = it->second;
      curr_key = key;
    }
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

  const std::string traits<ORBIT_PRINT_MODE>::name = "orbit_print_mode";

  const std::multimap<ORBIT_PRINT_MODE, std::vector<std::string> > traits<ORBIT_PRINT_MODE>::strval = {
    {ORBIT_PRINT_MODE::PROTO, {"PROTO", "Proto", "proto"} },
    {ORBIT_PRINT_MODE::FULL, {"FULL", "Full", "full"} }
  };

  ENUM_IO_DEF(ORBIT_PRINT_MODE)
  ENUM_JSON_IO_DEF(ORBIT_PRINT_MODE)

  jsonParser &to_json(const OrbitPrinterOptions &opt, jsonParser &json) {
    json.put_obj();
    json["indent_space"] = opt.indent_space;
    // just keep default delim
    json["prec"] = opt.prec;
    json[traits<COORD_TYPE>::name] = opt.coord_type;
    json[traits<ORBIT_PRINT_MODE>::name] = opt.orbit_print_mode;
    json["print_coordinates"] = opt.print_coordinates;
    json["print_equivalence_map"] = opt.print_equivalence_map;
    json["print_invariant_grp"] = opt.print_invariant_grp;
    json["sym_info_opt"] = opt.sym_info_opt;
    return json;
  }

  /// \brief Read from JSON
  void from_json(OrbitPrinterOptions &opt, const jsonParser &json) {
    json.get_if(opt.indent_space, "indent_space");
    // just keep default delim
    json.get_if(opt.prec, "prec");
    json.get_if(opt.coord_type, traits<COORD_TYPE>::name);
    json.get_if(opt.orbit_print_mode, traits<ORBIT_PRINT_MODE>::name);
    json.get_if(opt.print_coordinates, "print_coordinates");
    json.get_if(opt.print_equivalence_map, "print_equivalence_map");
    json.get_if(opt.print_invariant_grp, "print_invariant_grp");
    json.get_if(opt.sym_info_opt, "sym_info_opt");
  }

  OrbitPrinterOptions jsonConstructor<OrbitPrinterOptions>::from_json(const jsonParser &json) {
    OrbitPrinterOptions res;
    CASM::from_json(res, json);
    return res;
  }


  PrinterBase::PrinterBase(const OrbitPrinterOptions &_opt) :
    opt(_opt) {}

  void PrinterBase::coord_type(Log &out) {
    out << out.indent_str() << "COORD_MODE = " << opt.coord_type << std::endl << std::endl;
  }


  const std::string Printer<IntegralCluster>::element_name = "Clusters";

  void Printer<IntegralCluster>::print(const IntegralCluster &clust, Log &out) const {
    if(!out.print()) {
      return;
    }

    COORD_TYPE _mode = this->opt.coord_type;
    if(_mode == COORD_DEFAULT) {
      _mode = COORD_MODE::CHECK();
    }
    COORD_MODE printer_mode(_mode);
    if(_mode != INTEGRAL) {

      // calculate nice widths
      int prec = this->opt.prec;
      int width = prec;
      Eigen::Vector3d vec;
      out.ostream().precision(prec);
      out.ostream().flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      for(const auto &coord : clust) {
        if(_mode == CART) vec = coord.coordinate().cart();
        else if(_mode == FRAC) vec = coord.coordinate().frac();
        width = print_matrix_width(out, vec.transpose(), width);
      }

      // calculate nice widths
      Eigen::IOFormat format(prec, width + 1);
      for(const auto &coord : clust) {
        out << out.indent_str();
        coord.site().print(out, format);
        if(this->opt.delim) out << this->opt.delim;
        out << std::flush;
      }
    }
    else {
      // calculate nice widths
      int prec = 1;
      int width = prec;
      out.ostream().flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      for(const auto &coord : clust) {
        width = print_matrix_width(out, coord.unitcell().transpose(), width);
      }

      // print
      Eigen::IOFormat format(prec, width);
      for(const auto &coord : clust) {
        out << out.indent_str() << coord << " ";
        coord.site().occupant_dof().print(out);
        out << std::flush;
        if(this->opt.delim) out << this->opt.delim;
        out << std::flush;
      }
    }
  }

  ProtoFuncsPrinter::ProtoFuncsPrinter(ClexBasis const &_clex_basis, OrbitPrinterOptions const &_opt) :
    SitesPrinter(_opt),
    clex_basis(_clex_basis) {
    for(auto const &dofset : clex_basis.site_bases()) {
      for(BasisSet const &bset : dofset.second) {
        if(dofset.first != "occ"
           && bset.size()
           && bset[0]->type_name() != "Variable") {
          labelers.push_back(SubExpressionLabeler(bset.name(), "\\phi^{(" + dofset.first + ")}_%n_%l"));
        }
      }
    }
  }

  void print_site_basis_funcs(Structure const &prim,
                              ClexBasis const &clex_basis,
                              Log &out,
                              Index indent_space,
                              COORD_TYPE mode) {
    std::string indent(indent_space, ' ');

    std::ostream nullstream(0);
    COORD_MODE printer_mode(mode);
    std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > asym_unit;
    make_prim_periodic_asymmetric_unit(prim,
                                       CASM_TMP::ConstantFunctor<bool>(true),
                                       TOL,
                                       std::back_inserter(asym_unit),
                                       nullstream);


    for(auto const &dofset : clex_basis.site_bases()) {
      out << indent << indent << "Site basis functions for DoF \"" << dofset.first << "\":\n";
      for(Index no = 0; no < asym_unit.size(); no++) {
        out << indent << indent << "Asymmetric unit " << no + 1 << ":\n";
        for(Index ne = 0; ne < asym_unit[no].size(); ne++) {
          Index b = asym_unit[no][ne][0].sublat();
          out << indent << indent << "  Basis site " << b << ":\n"
              << "  ";
          if(printer_mode.check() == INTEGRAL) {
            out << indent << indent << asym_unit[no][ne][0] << ' ';
            asym_unit[no][ne][0].site().occupant_dof().print(out);
            out << std::flush;
          }
          else
            asym_unit[no][ne][0].site().print(out);

          out << "\n";
          if(dofset.second[b].size() == 0)
            out << "        [No site basis functions]\n\n";
          if(dofset.first == "occ") {
            for(Index f = 0; f < dofset.second[b].size(); f++) {

              BasisSet tbasis(dofset.second[b]);

              int s;
              std::vector<DoF::RemoteHandle> remote(1, DoF::RemoteHandle("occ", "s", prim.basis()[b].occupant_dof().ID()));
              remote[0] = s;
              tbasis.register_remotes(remote);

              for(s = 0; s < prim.basis()[b].occupant_dof().size(); s++) {
                if(s == 0)
                  out << "    ";
                out << "    \\phi_" << b << '_' << f << '[' << prim.basis()[b].occupant_dof()[s].name() << "] = "
                    << tbasis[f]->remote_eval();
                if(s + 1 == prim.basis()[b].occupant_dof().size())
                  out << "\n";
                else
                  out << ",   ";
              }
            }
          }
          else {
            //std::string formula;
            //bool relabel=false;
            for(Index f = 0; f < dofset.second[b].size(); f++) {
              if(dofset.second[b][f]) {
                out << "        ";
                if(dofset.second[b][f]->type_name() != "Variable") {
                  out << "\\phi^{" << (dofset.first) << "}_" << b << '_' << f << " = ";
                }
                out << dofset.second[b][f]->tex_formula() << "\n";
              }
            }
          }
        }
      }
    }
    out << "\n\n";


  }

  void write_site_basis_funcs(Structure const &prim,
                              ClexBasis const &clex_basis,
                              jsonParser &json) {


    //   "site_functions":[
    //     {
    //       "asym_unit": X,
    //       "sublat": 2,
    //       "basis": {
    //         "phi_b_0": {"Va":0.0, "O":1.0},
    //         "phi_b_1": {"Va":0.0, "O":1.0}
    //       }
    //     },
    //     ...
    //   ],

    jsonParser &sitef = json["site_functions"];
    sitef = jsonParser::array(prim.basis().size(), jsonParser::object());

    std::ostream nullstream(0);
    std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > asym_unit;
    make_prim_periodic_asymmetric_unit(prim,
                                       CASM_TMP::ConstantFunctor<bool>(true),
                                       TOL,
                                       std::back_inserter(asym_unit),
                                       nullstream);


    for(auto const &dofset : clex_basis.site_bases()) {
      for(Index no = 0; no < asym_unit.size(); no++) {

        for(Index ne = 0; ne < asym_unit[no].size(); ne++) {
          Index b = asym_unit[no][ne][0].sublat();
          sitef[b]["sublat"] = b;
          sitef[b]["asym_unit"] = no;

          if(dofset.second[b].size() == 0) {
            sitef[b][dofset.first]["basis"].put_null();
            continue;
          }


          for(Index f = 0; f < dofset.second[b].size(); f++) {
            std::stringstream fname;

            if(dofset.first == "occ") {
              fname << "\\phi_" << b << '_' << f;
              BasisSet tbasis(dofset.second[b]);
              int s;
              std::vector<DoF::RemoteHandle> remote(1, DoF::RemoteHandle("occ", "s", prim.basis()[b].occupant_dof().ID()));
              remote[0] = s;
              tbasis.register_remotes(remote);

              for(s = 0; s < prim.basis()[b].occupant_dof().size(); s++) {
                sitef[b][dofset.first]["basis"][fname.str()][prim.basis()[b].occupant_dof()[s].name()] = tbasis[f]->remote_eval();
              }
            }
            else {
              if(dofset.second[b][f]) {
                if(dofset.second[b][f]->type_name() != "Variable") {
                  fname << "\\phi^{" << (dofset.first) << "}_" << b << '_' << f;
                }
                else
                  fname << "var_" << b << '_' << f;
              }
              sitef[b][dofset.first]["basis"][fname.str()] = dofset.second[b][f]->tex_formula();
            }
          }
        }
      }
    }
  }

  // explicit template instantiations

#define PRINT_CLUST_INST(ITERATOR,INSERTER,PRINTER) \
  template void print_clust<ITERATOR,PRINTER>(ITERATOR begin, ITERATOR end, Log &out, PRINTER printer); \
  template jsonParser &write_clust<ITERATOR>(ITERATOR begin, ITERATOR end, jsonParser &json, PRINTER printer); \
  template jsonParser &write_clust<ITERATOR>(ITERATOR begin, ITERATOR end, jsonParser &json, PRINTER printer, const jsonParser &bspecs);

#define ORBIT_CONTAINER_INST(ITERATOR,INSERTER,ORBIT) \
  PRINT_CLUST_INST(ITERATOR,INSERTER,ProtoSitesPrinter) \
  PRINT_CLUST_INST(ITERATOR,INSERTER,FullSitesPrinter) \
  PRINT_CLUST_INST(ITERATOR,INSERTER,ProtoFuncsPrinter) \
  template void print_clust<ITERATOR>(ITERATOR begin, ITERATOR end, Log &out, const OrbitPrinterOptions &opt); \
  template INSERTER read_clust<INSERTER, typename ORBIT::SymCompareType>(\
    INSERTER result,\
    const jsonParser &json,\
    const Structure &prim,\
    const SymGroup& generating_grp,\
    const typename ORBIT::SymCompareType &sym_compare,\
    double xtal_tol);

#define _VECTOR_IT(ORBIT) std::vector<ORBIT>::iterator
#define _VECTOR_INSERTER(ORBIT) std::back_insert_iterator<std::vector<ORBIT> >

#define _SET_IT(ORBIT) std::set<ORBIT>::iterator
#define _SET_INSERTER(ORBIT) std::insert_iterator<std::set<ORBIT> >

#define ORBIT_VECTOR_INST(ORBIT) ORBIT_CONTAINER_INST(_VECTOR_IT(ORBIT),_VECTOR_INSERTER(ORBIT), ORBIT)
#define ORBIT_SET_INST(ORBIT) ORBIT_CONTAINER_INST(_SET_IT(ORBIT),_SET_INSERTER(ORBIT), ORBIT)

  ORBIT_VECTOR_INST(LocalIntegralClusterOrbit)
  ORBIT_VECTOR_INST(PrimPeriodicIntegralClusterOrbit)
  ORBIT_VECTOR_INST(ScelPeriodicIntegralClusterOrbit)
  ORBIT_VECTOR_INST(WithinScelIntegralClusterOrbit)

  ORBIT_SET_INST(LocalIntegralClusterOrbit)
  ORBIT_SET_INST(PrimPeriodicIntegralClusterOrbit)
  ORBIT_SET_INST(ScelPeriodicIntegralClusterOrbit)
  ORBIT_SET_INST(WithinScelIntegralClusterOrbit)


#define DIFFTRANS_VECTOR_INST(ORBIT) \
  PRINT_CLUST_INST(_VECTOR_IT(ORBIT), _VECTOR_INSERTER(ORBIT), PrototypePrinter<Kinetics::DiffusionTransformation>)

  DIFFTRANS_VECTOR_INST(PrimPeriodicDiffTransOrbit)
}
