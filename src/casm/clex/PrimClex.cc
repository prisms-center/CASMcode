#include "casm/clex/PrimClex.hh"

#include "casm/external/boost.hh"

#include "casm/misc/algorithm.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/ECIContainer.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clusterography/jsonClust.hh"
#include "casm/system/RuntimeLibrary.hh"
#include "casm/casm_io/SafeOfstream.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/app/AppIO.hh"
#include "casm/crystallography/Niggli.hh"

namespace CASM {
  //*******************************************************************************************
  //                                **** Constructors ****
  //*******************************************************************************************
  /// Initial construction of a PrimClex, from a primitive Structure
  PrimClex::PrimClex(const Structure &_prim, const Logging &logging) :
    Logging(logging),
    prim(_prim) {

    m_settings.set_crystallography_tol(TOL);

    _init();

    return;
  }


  //*******************************************************************************************
  /// Construct PrimClex from existing CASM project directory
  ///  - read PrimClex and directory structure to generate all its Supercells and Configurations, etc.
  PrimClex::PrimClex(const fs::path &_root, const Logging &logging):
    Logging(logging),
    root(_root),
    m_dir(_root),
    m_settings(_root),
    prim(read_prim(m_dir.prim())) {

    _init();

  }

  /// Initialization routines
  ///  - If !root.empty(), read all saved data to generate all Supercells and Configurations, etc.
  void PrimClex::_init() {

    log().construct("CASM Project");
    if(!root.empty()) {
      log() << "from: " << root << "\n" << std::endl;
    }

    std::vector<std::string> struc_mol_name = prim.get_struc_molecule_name();

    m_vacancy_allowed = false;
    for(int i = 0; i < struc_mol_name.size(); ++i) {
      if(is_vacancy(struc_mol_name[i])) {
        m_vacancy_allowed = true;
        m_vacancy_index = i;
      }
    }

    if(root.empty()) {
      return;
    }

    bool read_settings = false;
    bool read_composition = true;
    bool read_chem_ref = true;
    bool read_configs = true;
    bool clear_clex = false;

    refresh(read_settings, read_composition, read_chem_ref, read_configs, clear_clex);
  }

  /// \brief Reload PrimClex data from settings
  ///
  /// \param read_settings Read project_settings.json and plugins
  /// \param read_composition Read composition_axes.json
  /// \param read_chem_ref Read chemical_reference.json
  /// \param read_configs Read SCEL and config_list.json
  /// \param clear_clex Clear stored orbitrees, clexulators, and eci
  ///
  /// - This does not check if what you request will cause problems.
  /// - ToDo: refactor into separate functions
  ///
  void PrimClex::refresh(bool read_settings,
                         bool read_composition,
                         bool read_chem_ref,
                         bool read_configs,
                         bool clear_clex) {

    log().custom("Load project data");

    if(read_settings) {
      try {
        m_settings = ProjectSettings(root, *this);
      }
      catch(std::exception &e) {
        err_log().error("reading project_settings.json");
        err_log() << "file: " << m_dir.project_settings() << "\n" << std::endl;
      }
    }

    if(read_composition) {
      m_has_composition_axes = false;
      auto comp_axes = m_dir.composition_axes();

      try {
        if(fs::is_regular_file(comp_axes)) {
          log() << "read: " << comp_axes << "\n";

          CompositionAxes opt(comp_axes);

          if(opt.has_current_axes) {
            m_has_composition_axes = true;
            m_comp_converter = opt.curr;
          }
        }
      }
      catch(std::exception &e) {
        err_log().error("reading composition_axes.json");
        err_log() << "file: " << comp_axes << "\n" << std::endl;
      }
    }

    if(read_chem_ref) {

      // read chemical reference
      m_chem_ref.reset();
      auto chem_ref_path = m_dir.chemical_reference(m_settings.default_clex().calctype, m_settings.default_clex().ref);

      try {
        if(fs::is_regular_file(chem_ref_path)) {
          log() << "read: " << chem_ref_path << "\n";
          m_chem_ref = notstd::make_cloneable<ChemicalReference>(read_chemical_reference(chem_ref_path, prim, lin_alg_tol()));
        }
      }
      catch(std::exception &e) {
        err_log().error("reading chemical_reference.json");
        err_log() << "file: " << chem_ref_path << "\n" << std::endl;
      }
    }

    if(read_configs) {

      supercell_list.clear();

      try {
        // read supercells
        if(fs::is_regular_file(m_dir.SCEL())) {
          log() << "read: " << m_dir.SCEL() << "\n";
          fs::ifstream scel(m_dir.SCEL());
          read_supercells(scel);
        }
      }
      catch(std::exception &e) {
        err_log().error("reading SCEL");
        err_log() << "file: " << m_dir.SCEL() << "\n" << std::endl;
      }

      try {
        // read config_list
        if(fs::is_regular_file(get_config_list_path())) {
          log() << "read: " << get_config_list_path() << "\n";
          read_config_list();
        }
      }
      catch(std::exception &e) {
        err_log().error("reading config_list.json");
        err_log() << "file: " << m_dir.config_list() << "\n" << std::endl;
      }
    }

    if(clear_clex) {
      m_nlist.reset();
      m_orbitree.clear();
      m_clexulator.clear();
      m_eci.clear();
      log() << "refresh cluster expansions\n";
    }

    log() << std::endl;
  }

  //*******************************************************************************************
  // **** Accessors ****
  //*******************************************************************************************

  /// Return project name
  std::string PrimClex::name() const {
    return settings().name();
  }

  //*******************************************************************************************
  // ** Directory path accessors **
  //*******************************************************************************************
  /// Return casm project directory path
  fs::path PrimClex::get_path() const {
    return root;
  }

  //*******************************************************************************************
  /// Return supercell directory path
  fs::path PrimClex::get_path(const Index &scel_index) const {
    return root / "training_data" / supercell_list[scel_index].get_name();
  }

  //*******************************************************************************************
  /// Return configuration directory path
  fs::path PrimClex::get_path(const Index &scel_index, const Index &config_index) const {
    return get_path(scel_index) / supercell_list[scel_index].get_config(config_index).get_id();
  }

  //*******************************************************************************************
  /// Return config_list.json file path
  fs::path PrimClex::get_config_list_path() const {
    return root / ".casm" / "config_list.json";
  }


  // ** Current settings accessors **

  // ** Composition accessors **

  //*******************************************************************************************
  /// const Access CompositionConverter object
  bool PrimClex::has_composition_axes() const {
    return m_has_composition_axes;
  }

  //*******************************************************************************************
  /// const Access CompositionConverter object
  const CompositionConverter &PrimClex::composition_axes() const {
    return m_comp_converter;
  }

  // ** Chemical reference **

  //*******************************************************************************************
  /// check if ChemicalReference object initialized
  bool PrimClex::has_chemical_reference() const {
    return static_cast<bool>(m_chem_ref);
  }

  //*******************************************************************************************
  /// const Access ChemicalReference object
  const ChemicalReference &PrimClex::chemical_reference() const {
    return *m_chem_ref;
  }


  // ** Prim and Orbitree accessors **

  //*******************************************************************************************
  /// const Access to primitive Structure
  const Structure &PrimClex::get_prim() const {
    return prim;
  }

  //*******************************************************************************************

  PrimNeighborList &PrimClex::nlist() const {

    // lazy neighbor list generation
    if(!m_nlist) {

      // construct nlist
      m_nlist = notstd::make_cloneable<PrimNeighborList>(
                  settings().nlist_weight_matrix(),
                  settings().nlist_sublat_indices().begin(),
                  settings().nlist_sublat_indices().end()
                );
    }

    return *m_nlist;
  }

  //*******************************************************************************************
  /// returns true if vacancy are an allowed species
  bool PrimClex::vacancy_allowed() const {
    return m_vacancy_allowed;
  }

  //*******************************************************************************************
  /// returns the index of vacancies in composition vectors
  Index PrimClex::vacancy_index() const {
    return m_vacancy_index;
  }


  // ** Supercell and Configuration accessors **

  //*******************************************************************************************
  /// const Access entire supercell_list
  boost::container::stable_vector<Supercell> &PrimClex::get_supercell_list() {
    return supercell_list;
  };

  //*******************************************************************************************
  /// const Access entire supercell_list
  const boost::container::stable_vector<Supercell> &PrimClex::get_supercell_list() const {
    return supercell_list;
  };

  //*******************************************************************************************
  /// const Access supercell by index
  const Supercell &PrimClex::get_supercell(Index i) const {
    return supercell_list[i];
  };

  //*******************************************************************************************
  /// Access supercell by index
  Supercell &PrimClex::get_supercell(Index i) {
    return supercell_list[i];
  };

  //*******************************************************************************************
  /// const Access supercell by name
  const Supercell &PrimClex::get_supercell(std::string scellname) const {
    Index index;
    if(!contains_supercell(scellname, index)) {
      err_log().error("Accessing supercell");
      err_log() << "supercell '" << scellname << "' not found." << std::endl;
      throw std::invalid_argument("Error in PrimClex::get_supercell(std::string scellname) const: Not found");
    }
    return supercell_list[index];
  };

  //*******************************************************************************************
  /// Access supercell by name
  Supercell &PrimClex::get_supercell(std::string scellname) {
    return const_cast<Supercell &>(static_cast<const PrimClex &>(*this).get_supercell(scellname));
  }

  //*******************************************************************************************
  /// Access supercell by Lattice, adding if necessary
  Supercell &PrimClex::get_supercell(const Lattice &lat) {
    return get_supercell(add_supercell(lat));
  }

  //*******************************************************************************************
  /// access configuration by name (of the form "scellname/[NUMBER]", e.g., ("SCEL1_1_1_1_0_0_0/0")
  const Configuration &PrimClex::configuration(const std::string &configname) const {
    auto res = Configuration::split_name(configname);

    try {
      return get_supercell(res.first).get_config_list().at(res.second);
    }
    catch(std::out_of_range &e) {
      err_log().error("Invalid config index");
      err_log() << "ERROR: In PrimClex::configuration(), configuration index out of range\n";
      err_log() << "configname: " << configname << "\n";
      err_log() << "index: " << res.second << "\n";
      err_log() << "config_list.size(): " << get_supercell(res.first).get_config_list().size() << "\n";
      throw e;
    }
  }

  //*******************************************************************************************

  Configuration &PrimClex::configuration(const std::string &configname) {
    return const_cast<Configuration &>(static_cast<const PrimClex &>(*this).configuration(configname));
  }

  //*******************************************************************************************
  /// Configuration iterator: begin
  PrimClex::config_iterator PrimClex::config_begin() {
    if(supercell_list.size() == 0 || supercell_list[0].get_config_list().size() > 0)
      return config_iterator(this, 0, 0);
    return ++config_iterator(this, 0, 0);
  }

  //*******************************************************************************************
  /// Configuration iterator: end
  PrimClex::config_iterator PrimClex::config_end() {
    return config_iterator(this, supercell_list.size(), 0);
  }

  //*******************************************************************************************
  /// const Configuration iterator: begin
  PrimClex::config_const_iterator PrimClex::config_begin() const {
    if(supercell_list.size() == 0 || supercell_list[0].get_config_list().size() > 0)
      return config_const_iterator(this, 0, 0);
    return ++config_const_iterator(this, 0, 0);
  }

  //*******************************************************************************************
  /// const Configuration iterator: end
  PrimClex::config_const_iterator PrimClex::config_end() const {
    return config_const_iterator(this, supercell_list.size(), 0);
  }

  //*******************************************************************************************
  /// const Configuration iterator: begin
  PrimClex::config_const_iterator PrimClex::config_cbegin() const {
    if(supercell_list.size() == 0 || supercell_list[0].get_config_list().size() > 0)
      return config_const_iterator(this, 0, 0);
    return ++config_const_iterator(this, 0, 0);
  }

  //*******************************************************************************************
  /// const Configuration iterator: end
  PrimClex::config_const_iterator PrimClex::config_cend() const {
    return config_const_iterator(this, supercell_list.size(), 0);
  }

  //*******************************************************************************************
  /// Configuration iterator: begin
  PrimClex::config_iterator PrimClex::selected_config_begin() {
    //std::cout << "BEGINNING SELECTED CONFIG ITERATOR\n"
    //          << "supercell_list.size() is " << supercell_list.size() << "\n";

    if(supercell_list.size() == 0 || (supercell_list[0].get_config_list().size() > 0 && supercell_list[0].get_config(0).selected()))
      return config_iterator(this, 0, 0, true);
    return ++config_iterator(this, 0, 0, true);
  }

  //*******************************************************************************************
  /// Configuration iterator: end
  PrimClex::config_iterator PrimClex::selected_config_end() {
    return config_iterator(this, supercell_list.size(), 0, true);
  }

  //*******************************************************************************************
  /// const Configuration iterator: begin
  PrimClex::config_const_iterator PrimClex::selected_config_cbegin() const {
    if(supercell_list.size() == 0 || (supercell_list[0].get_config_list().size() > 0 && supercell_list[0].get_config(0).selected()))
      return config_const_iterator(this, 0, 0, true);
    return ++config_const_iterator(this, 0, 0, true);
  }

  //*******************************************************************************************
  /// const Configuration iterator: end
  PrimClex::config_const_iterator PrimClex::selected_config_cend() const {
    return config_const_iterator(this, supercell_list.size(), 0, true);
  }


  //*******************************************************************************************
  // **** IO ****
  //*******************************************************************************************
  /**
   * Re-write config_list.json, excluding specified scel
   */

  void PrimClex::write_config_list(std::set<std::string> scel_to_delete) {

    if(supercell_list.size() == 0) {
      fs::remove(get_config_list_path());
      return;
    }

    jsonParser json;

    if(fs::exists(get_config_list_path())) {
      json.read(get_config_list_path());
    }
    else {
      json.put_obj();
    }

    for(Index s = 0; s < supercell_list.size(); s++) {
      if(scel_to_delete.count(supercell_list[s].get_name())) {
        json["supercells"].erase(supercell_list[s].get_name());
      }
      else {
        supercell_list[s].write_config_list(json);
      }
    }

    SafeOfstream file;
    file.open(get_config_list_path());
    json.print(file.ofstream());
    file.close();

    return;
  }

  // **** Operators ****


  // **** Unorganized mess of functions ... ****

  /// \brief Generate supercells of a certain volume and shape and store them in the array of supercells
  ///
  /// \param enum_props An ScelEnumProps instance, see constructor for details
  ///
  void PrimClex::generate_supercells(const ScelEnumProps &enum_props) {
    ScelEnumByProps e(*this, enum_props);
    for(auto it = e.begin(); it != e.end(); ++it) {}
    return;
  }

  //*******************************************************************************************
  /**
   *  add a supercell if it doesn't already exist
   *    return the index of the supercell in supercell_list
   *
   *    TODO:
   *    Check to see if superlat is linear combination of previously
   *    enumerated Supercell. You'll have to transform the coordinates
   *    if you're trying to import something with a supercell that
   *    already exists. Should return operations involved.
   *    If supercell doesn't exist add it as it is without transforming
   */
  //*******************************************************************************************
  Index PrimClex::add_canonical_supercell(const Lattice &superlat) {
    if(!superlat.is_supercell_of(prim.lattice())) {
      std::cerr << "ERROR in PrimClex::add_canonical_supercell()." << std::endl
                << "  Passed a Supercell Lattice that is not a superlattice of PRIM lattice\n" << std::endl;
      assert(0);
      exit(1);
    }

    // temporary version, just check transf_mat equivalence
    // Does this check for equivalent supercells with different transformation matrices? Seems like it should
    // Insert second loop that goes over a symmetry operation list and applies it to the transformation matrix
    Supercell scel(this, superlat);
    scel.set_id(supercell_list.size());
    for(Index i = 0; i < supercell_list.size(); i++) {
      if(supercell_list[i].get_transf_mat() == scel.get_transf_mat())
        return i;
    }

    // if not already existing, add it
    supercell_list.push_back(scel);
    return supercell_list.size() - 1;
  }
  //*******************************************************************************************
  /**
   *  add a supercell if it doesn't already exist
   *  return the index of the supercell in supercell_list
   *
   *  Unlike the initial version above this, this routine
   *  will use an Array of SymOps (probably the point group of the lattice)
   *  to determine whether the provided superlat is equivalent to a Supercell in the
   *  PrimClex list. It will then populate the transformation matrix that
   *  maps the passed lattice onto the one on the list (i.e. the one on on
   *  the list is treated as 'primitive'). If it doesn't map, then the
   *  matrix gets flattened to zeros.
   *
   *  This routine should replace the original one once it converges
   *  to something people can agree on.
   */
  //*******************************************************************************************
  Index PrimClex::add_supercell(const Lattice &superlat) {

    return add_canonical_supercell(canonical_equivalent_lattice(superlat, prim.point_group(), crystallography_tol()));

  }

  //*******************************************************************************************
  void PrimClex::print_supercells(std::set<std::string> scel_to_delete) const {


    // also write supercells/supercell_path directories with LAT files
    try {
      fs::create_directory(get_path() / "training_data");
    }
    catch(const fs::filesystem_error &ex) {
      std::cerr << "Error in PrimClex::print_supercells()." << std::endl;
      std::cerr << ex.what() << std::endl;
    }

    fs::ofstream scelfile(get_path() / "training_data" / "SCEL");

    print_supercells(scelfile, scel_to_delete);

    for(Index i = 0; i < supercell_list.size(); i++) {
      try {
        if(!scel_to_delete.count(supercell_list[i].get_name())) {
          fs::create_directory(supercell_list[i].get_path());

          fs::path latpath = supercell_list[i].get_path() / "LAT";
          if(!fs::exists(latpath)) {
            fs::ofstream latfile(latpath);
            supercell_list[i].get_real_super_lattice().print(latfile);
          }
        }
      }
      catch(const fs::filesystem_error &ex) {
        std::cerr << "Error in PrimClex::print_supercells()." << std::endl;
        std::cerr << ex.what() << std::endl;
      }
    }
  }

  //*******************************************************************************************
  void PrimClex::print_supercells(std::ostream &stream, std::set<std::string> scel_to_delete) const {
    for(Index i = 0; i < supercell_list.size(); i++) {
      if(!scel_to_delete.count(supercell_list[i].get_name())) {
        stream << "Supercell Name: '" << supercell_list[i].get_name() << "' Number: " << i << " Volume: " << supercell_list[i].get_transf_mat().determinant() << "\n";
        stream << "Supercell Transformation Matrix: \n";
        stream << supercell_list[i].get_transf_mat();
        stream << "\n";
      }
    }
  }

  //*******************************************************************************************
  void PrimClex::read_supercells(std::istream &stream) {
    // expect a file with format:
    //
    // Supercell Number: 0 Volume: 1
    // Supercell Transformation Matrix:
    //  1 0 0
    //  0 1 0
    //  0 0 1
    //
    // Supercell Number: 1 Volume: 2
    // Supercell Transformation Matrix:
    //  1 0 -1
    //  0 1 0
    //  0 0 2

    Eigen::Matrix3d mat;

    // check for non-canonical
    std::vector<std::pair<Lattice, Lattice> > non_canon;

    // read && sort
    std::set<Supercell> in_scel;

    std::string s;
    while(!stream.eof()) {
      std::getline(stream, s);
      if(s[0] == 'S') {
        std::getline(stream, s);
        stream >> mat;

        Lattice lat(prim.lattice().lat_column_mat()*mat);
        in_scel.emplace(this, lat);

        Lattice canon = canonical_equivalent_lattice(
                          lat,
                          get_prim().point_group(),
                          crystallography_tol());

        if(!almost_equal(canon.lat_column_mat(), lat.lat_column_mat(), crystallography_tol())) {
          non_canon.push_back(std::make_pair(lat, canon));
        }
      }
    }

    // insert sorted
    for(const auto &scel : in_scel) {
      add_canonical_supercell(scel.get_real_super_lattice());
    }

    if(non_canon.size()) {
      log() << std::endl;

      log().warning("Non-canonical supercells");
      log() << "A bug in a previous version of CASM resulted in the generation \n"
            "of the following non-canonical supercells: \n";
      for(const auto &val : non_canon) {
        log() << "  existing name: " << generate_name(calc_transf_mat(val.first))
              << "  canonical name: " << generate_name(calc_transf_mat(val.second)) << "\n";
      }
      log() << std::endl;
      log() << "To prevent errors, please use 'casm import' to convert existing \n"
            "configurations and data.\n";
      log() << std::endl;

    }
  }

  //*******************************************************************************************
  /**
   *   Read the config_list file at 'file_name' and adds config
   *   to supercell, assuming it is already canonical.
   *   If 'print_dirs', call Supercell::print_clex_configuration()
   *   for each config to be made.
   */
  //*******************************************************************************************
  void PrimClex::read_config_list() {

    jsonParser json(get_config_list_path());

    if(!json.contains("supercells")) {
      return;
    }

    for(Index i = 0; i < supercell_list.size(); i++) {
      supercell_list[i].read_config_list(json);
    }
  }

  //*******************************************************************************************
  /*
   * Run through all the supercells and add up how many configurations are selected
   * in each one of them.
   */
  //*******************************************************************************************

  int PrimClex::amount_selected() const {
    int amount_selected = 0;
    for(Index s = 0; s < supercell_list.size(); s++) {
      amount_selected += supercell_list[s].amount_selected();
    }
    return amount_selected;
  }

  //*******************************************************************************************
  bool PrimClex::contains_supercell(std::string scellname, Index &index) const {
    for(Index i = 0; i < supercell_list.size(); i++) {
      if(supercell_list[i].get_name() == scellname) {
        index = i;
        return true;
      }
    }
    index = supercell_list.size();
    return false;

  };

  bool PrimClex::contains_supercell(const Supercell &scel) const {
    Index tmp;
    return contains_supercell(scel, tmp);
  }

  bool PrimClex::contains_supercell(const Supercell &scel, Index &index) const {
    return contains_supercell(scel.get_name(), index);
  }

  //*******************************************************************************************
  Eigen::Matrix3i PrimClex::calc_transf_mat(const Lattice &superlat) const {
    Eigen::Matrix3d ttrans = prim.lattice().inv_lat_column_mat() * superlat.lat_column_mat();
    if(!is_integer(ttrans, TOL)) {
      std::cerr << "Error in PrimClex::calc_transf_mat(const Lattice &superlat)" << std::endl
                << "  Bad supercell, the transformation matrix is not integer. Exiting!" << std::endl;
      exit(1);
    }
    return iround(ttrans);
  }

  //*******************************************************************************************

  /// \brief Sets the composition axes
  ///
  /// Also:
  /// - updates all configuration references,
  /// - writes the updated configuration info
  /// - does not update composition_axes file
  void PrimClex::set_composition_axes(const CompositionConverter &_converter) {

    m_comp_converter = _converter;
    m_has_composition_axes = true;

  }

  //*******************************************************************************************
  /*
    /// Delete 'properties.ref_state.X.json' files,
    /// Then call 'clear_reference_properties'
    void PrimClex::clear_reference_states() {
      for(int i = 0; i < composition_axes().independent_compositions() + 1; i++) {
        fs::remove(dir().ref_state(settings().calctype(), settings().ref(), i));
      }
      generate_references();
    }
  */
  //*******************************************************************************************
  /*
    /// Sets the root reference state to be the calculated properties of the chosen config
    /// Does not call 'generate_references' so written files will be out-of-date until you do so!!!
    void PrimClex::set_reference_state(int refid, const Configuration &config) {

      //std::cout << "begin set_reference_state()" << std::endl;

      std::string reason_invalid;
      if(!valid_reference_state(refid, config, reason_invalid)) {
        std::cerr << "Error in PrimClex::set_reference_state." << std::endl
                  << "  Could not use  SCEL: " << config.get_supercell().get_name() << " ConfigID: " << config.get_id() << "  for reference state: " << refid << std::endl
                  << "  " << reason_invalid << std::endl;
        exit(1);
      }

      jsonParser json;

      json["supercell_name"] = config.get_supercell().get_name();
      json["configid"] = config.get_id();

      json["param_composition"] = config.get_param_composition();
      json["ref_state"] = config.calc_properties();

      //std::cout << "Set refid: " << refid << std::endl;
      //std::cout << "Reference State:\n" << json << "\n" << std::endl;

      json.write(dir().ref_state(settings().calctype(), settings().ref(), refid));

      //std::cout << "finish set_reference_state()" << std::endl;
    }
  */
  //*******************************************************************************************
  /*
    /// Check that it is valid to use 'config' as reference state 'refid', returns bool and if false, sets 'reason_invalid'
    ///   Currently checks:
    ///     1) that the necessary properties have been calculated,
    ///     2) that the same Configuration is not being used twice
    ///   Needs to check that reference states span composition space
    bool PrimClex::valid_reference_state(int refid, const Configuration &config, std::string &reason_invalid) const {

      // Check that the Configuration has all the curr_property calculated
      for(Index i = 0; i < get_curr_property().size(); i++) {
        if(!config.calc_properties().contains(get_curr_property()[i])) {
          reason_invalid = "You are attempting to use a Configuration for which the property '" + get_curr_property()[i] + "' has not been calculated.";
          return false;
        }
      }

      // Check that a Configuration is not being used for multiple reference states
      for(int i = 0; i < composition_axes().independent_compositions() + 1; i++) {
        if(i == refid)
          continue;

        if(!fs::exists(dir().ref_state(settings().calctype(), settings().ref(), i)))
          continue;

        jsonParser json(dir().ref_state(settings().calctype(), settings().ref(), i));

        if(json["configid"] == "custom" || json["supercell_name"] == "custom")
          continue;

        if(json["configid"] == config.get_id() && json["supercell_name"] == config.get_supercell().get_name()) {
          reason_invalid =  "You are attempting to use the same Configuration for >1 reference state and I can't allow that.";
          return false;
        }
      }

      // Here I should check that references span the necessary composition space, but I'm not yet

      //   First read all reference states that have already been set, then check that this new one is not a linear combination of those

      return true;
    }
  */
  //*******************************************************************************************
  /*
    /// find calculated configurations closest to
    /// [0, 0, 0, ...], [1, 0, 0, ...], [0, 1, 0, ...], [0, 0, 1, ...], ...
    /// and set them as the root reference states, also calls generate_references
    /// Clears refrence states and properties whether or not it succeeds
    void PrimClex::set_reference_state_auto() {

      //std::cout << "begin set_reference_state_auto()" << std::endl;

      /// find calculated configurations closest to
      /// [0, 0, 0, ...], [1, 0, 0, ...], [0, 1, 0, ...], [0, 0, 1, ...], ...

      // Clear current references
      clear_reference_states();

      int Naxes = composition_axes().independent_compositions();

      //std::cout << "Naxes: " << Naxes << std::endl;

      Eigen::VectorXd target = Eigen::VectorXd::Zero(Naxes);

      //std::cout << "Ref: " << 0 << "  target: " << target.transpose() << std::endl;
      set_reference_state(0, closest_calculated_config(target));

      for(int i = 0; i < Naxes; i++) {

        target(i) = 1.0;
        //std::cout << "Ref: " << i + 1 << "  target: " << target.transpose() << std::endl;
        set_reference_state(i + 1, closest_calculated_config(target));
        target(i) = 0.0;
      }

      //std::cout << "generate references" << std::endl;
      generate_references();

      //std::cout << "finish set_reference_state_auto()" << std::endl;

    }
  */

  //*******************************************************************************************
  /*
    /// Clear 'reference' and 'delta' properties from all Configurations
    /// Re-write all Configurations, updating:
    ///   param_composition.json
    ///   properties.calc.json
    ///   properties.ref.json
    ///   properties.delta.json
    void PrimClex::generate_references() {

      //std::cout << "begin PrimClex::generate_references()" << std::endl;

      for(config_iterator it = config_begin(); it != config_end(); ++it) {
        it->generate_reference();
      }

      write_config_list();

      //std::cout << "finish PrimClex::generate_references()" << std::endl;

    }
  */

  //*******************************************************************************************
  /// private:

  //*******************************************************************************************
  /*
    /// Return the configuration closest in param_composition to the target_param_comp
    ///   Tie break returns configuration in smallest supercell (first found at that size)
    const Configuration &PrimClex::closest_calculated_config(const Eigen::VectorXd &target_param_comp) const {

      //std::cout << "begin closest_calculated_config()" << std::endl;

      /// return reference to Configuration with param_comp closest to target_param_comp
      ///   tie break goes to first Configuration with fewest atoms
      ///
      ///   must be Configurations for which the curr_properties have been calculated

      Eigen::VectorXd param_comp;
      Eigen::VectorXd closest_comp;

      double curr_dist;
      double close_dist = -1;
      Index close_super = -1;
      Index close_config;
      Index close_size;

      for(Index i = 0; i < supercell_list.size(); i++) {
        for(Index j = 0; j < supercell_list[i].get_config_list().size(); j++) {

          const Configuration &config = supercell_list[i].get_config(j);

          // check if the config has been calculated
          //std::cout << "\n\nScell: " << supercell_list[i].get_name() << "  Config: " << j << "\n" << config.calc_properties() << std::endl;

          if(config.calc_properties().size() == 0)
            continue;

          curr_dist = (target_param_comp - config.get_param_composition()).norm();

          if(!valid_index(close_super) ||
             (almost_equal(curr_dist, close_dist, TOL) && config.size() < close_size) ||
             (curr_dist < close_dist)) {
            close_super = i;
            close_config = j;
            close_dist = curr_dist;
            close_size = config.size();
          }

        } // for j
      } // for i

      if(!valid_index(close_super)) {
        std::cerr << "Error in PrimClex::closest_calculated_config" << std::endl
                  << "  Could not find a calculated Configuration." << std::endl;
        exit(1);
      }

      //std::cout << "Closest Calculated:  SCEL: " << supercell_list[close_super].get_name() << "  Config: " <<
      //          close_config << "  ParamComp: " << supercell_list[close_super].get_config(close_config).get_param_composition().transpose() << std::endl;

      //std::cout << "finish closest_calculated_config()" << std::endl;

      return supercell_list[close_super].get_config(close_config);
    }
  */

  //*******************************************************************************************

  Eigen::MatrixXd PrimClex::shift_vectors() const {
    Eigen::MatrixXd tau(prim.basis.size(), 3);
    for(int i = 0; i < prim.basis.size(); i++) {
      tau.row(i) = prim.basis[i].const_cart().transpose();
    }
    return tau;
  }

  //*******************************************************************************************
  /// const Access to global orbitree
  bool PrimClex::has_orbitree(const ClexDescription &key) const {
    auto it = m_orbitree.find(key);
    if(it == m_orbitree.end()) {
      if(!fs::exists(dir().clust(key.bset))) {
        return false;
      }
    }
    return true;

  };

  //*******************************************************************************************

  const SiteOrbitree &PrimClex::orbitree(const ClexDescription &key) const {

    auto it = m_orbitree.find(key);
    if(it == m_orbitree.end()) {
      it = m_orbitree.insert(std::make_pair(key, SiteOrbitree(get_prim().lattice(), crystallography_tol()))).first;
      SiteOrbitree &tree = it->second;

      // these could be specified via settings
      tree.min_num_components = 2;
      tree.min_length = 0.0001;

      from_json(jsonHelper(tree, prim), jsonParser(dir().clust(key.bset)));
      tree.generate_clust_bases();

    }

    return it->second;

  }

  //*******************************************************************************************

  bool PrimClex::has_clexulator(const ClexDescription &key) const {
    auto it = m_clexulator.find(key);
    if(it == m_clexulator.end()) {
      if(!fs::exists(dir().clexulator_src(settings().name(), key.bset))) {
        return false;
      }
    }
    return true;
  }

  //*******************************************************************************************

  Clexulator PrimClex::clexulator(const ClexDescription &key) const {

    auto it = m_clexulator.find(key);
    if(it == m_clexulator.end()) {

      if(!fs::exists(dir().clexulator_src(settings().name(), key.bset))) {
        throw std::runtime_error(
          std::string("Error loading clexulator ") + key.bset + ". No basis functions exist.");
      }

      try {
        it = m_clexulator.insert(
               std::make_pair(key, Clexulator(settings().name() + "_Clexulator",
                                              dir().clexulator_dir(key.bset),
                                              nlist(),
                                              log(),
                                              settings().compile_options(),
                                              settings().so_options()))).first;
      }
      catch(std::exception &e) {
        // not sure why this fails...
        // log() << "Error constructing Clexulator. Current settings: \n" << std::endl;
        // settings().print_compiler_settings_summary(log());

        std::cout << "Error constructing Clexulator. Current settings: \n" << std::endl;
        Log tlog(std::cout);
        settings().print_compiler_settings_summary(tlog);
        throw;
      }
    }
    return it->second;
  }

  //*******************************************************************************************

  bool PrimClex::has_eci(const ClexDescription &key) const {

    auto it = m_eci.find(key);
    if(it == m_eci.end()) {
      return fs::exists(dir().eci(key.property, key.calctype, key.ref, key.bset, key.eci));
    }
    return true;
  }

  //*******************************************************************************************

  const ECIContainer &PrimClex::eci(const ClexDescription &key) const {

    auto it = m_eci.find(key);
    if(it == m_eci.end()) {
      fs::path eci_path = dir().eci(key.property, key.calctype, key.ref, key.bset, key.eci);
      if(!fs::exists(eci_path)) {
        throw std::runtime_error(
          std::string("Error loading ECI. eci.json does not exist.\n")
          + "  Expected at: " + eci_path.string());
      }

      it = m_eci.insert(std::make_pair(key, read_eci(eci_path))).first;
    }
    return it->second;
  }

  //*******************************************************************************************
  /// \brief Make orbitree. For now specifically global.
  ///
  /// \param prim Primitive Structure. non-const due to Structure::set_site_internals.
  /// \param json bspecs.json file
  SiteOrbitree make_orbitree(Structure &prim, const jsonParser &json, double _tol) {

    try {

      SiteOrbitree tree(prim.lattice(), _tol);
      tree.set_bspecs(json);
      // --- first generate global orbitree -------------


      //global_orbitree.read_CSPECS(in_clust);
      tree.min_num_components = 2;
      tree.min_length = -1.5 * _tol;

      tree.max_length.clear();
      auto update_max_length = [&](int branch, double max_length) {
        while(branch > tree.max_length.size() - 1) {
          tree.max_length.push_back(0.0);
        }
        tree.max_length[branch] = max_length;
      };

      for(auto it = json["orbit_branch_specs"].cbegin(); it != json["orbit_branch_specs"].cend(); ++it) {
        update_max_length(std::stoi(it.name()), it->find("max_length")->get<double>());
      }
      tree.max_num_sites = tree.max_length.size() - 1;

      tree.generate_orbitree(prim);

      // --- then add custom orbits --------------------

      if(json.contains("orbit_specs")) {
        bool verbose = false;
        tree.read_custom_clusters_from_json(json["orbit_specs"], prim, prim.factor_group(), verbose);
      }
      tree.collect_basis_info(prim);

      return tree;
    }
    catch(...) {
      std::cerr << "Error in make_orbitree, with JSON input: \n";
      std::cerr << json << "\n";
      throw;
    }

  }

  void set_nlist_ind(const Structure &prim, SiteOrbitree &tree, const PrimNeighborList &nlist, double xtal_tol) {

    //For each site we encounter we access the appropriate slot in the neighbor list and append all other sites
    //to it in the form of UnitCellCoords

    Array<Index> clust_nlist_inds;

    const auto &sublat_indices = nlist.sublat_indices();

    Index N_sublat = sublat_indices.size();

    //branches
    for(Index i = 0; i < tree.size(); i++) {
      //orbits
      for(Index j = 0; j < tree[i].size(); j++) {
        //clusters
        for(Index k = 0; k < tree[i][j].size(); k++) {

          clust_nlist_inds.resize(tree[i][j][k].size());

          //sites
          for(Index l = 0; l < tree[i][j][k].size(); l++) {

            //tuccl corresponds to a particular site we're looking at
            UnitCellCoord tuccl(tree[i][j][k][l], prim, xtal_tol);

            //neighbor sites
            for(Index b = 0; b < tree[i][j][k].size(); b++) {

              //tuccb corresponds to a site that neighbors tuccl
              UnitCellCoord tuccb(tree[i][j][k][b], prim, xtal_tol);
              UnitCell delta = tuccb.unitcell() - tuccl.unitcell();

              auto unitcell_index = find_index(nlist, delta);
              if(unitcell_index == nlist.size()) {
                std::cerr << "Error generating unitcell index." << std::endl;
                std::cerr << "  Did not find unitcell: " << delta.transpose() << " in the prim nlist." << std::endl;
                exit(1);
              }

              auto sublat_index = find_index(sublat_indices, tuccb.sublat());
              if(sublat_index == sublat_indices.size()) {
                std::cerr << "Error generating sublat index" << std::endl;
                std::cerr << "  Did not find sublat: " << tuccb.sublat() << " in the nlist sublattice indices: " << jsonParser(sublat_indices) << std::endl;
                exit(1);
              }

              clust_nlist_inds[b] = unitcell_index * N_sublat + sublat_index;

              // //If delta is not already in nlist, add it (the value of clust_nlist_inds[b] makes sense after the push_back)
              // if(clust_nlist_inds[b] == nlist.size()) {
              //  nlist.push_back(delta);
              // }
            }

            // we set the first set of nlist_inds as the current indices
            if(l == 0) {
              tree[i][j][k].set_nlist_inds(clust_nlist_inds);
            }

            // save any other unique ones that we find
            Index t(0);
            for(t = 0; t < tree[i][j][k].trans_nlists().size(); t++) {
              if(clust_nlist_inds.all_in(tree[i][j][k].trans_nlist(t)))
                break;
            }
            if(t == tree[i][j][k].trans_nlists().size())
              tree[i][j][k].add_trans_nlist(clust_nlist_inds);
          }
        }
      }
    }
  }

  //*******************************************************************************************
  /// \brief Print clexulator
  void print_clexulator(const Structure &prim,
                        SiteOrbitree &tree,
                        const PrimNeighborList &nlist,
                        std::string class_name,
                        std::ostream &stream,
                        double xtal_tol) {

    set_nlist_ind(prim, tree, nlist, xtal_tol);

    DoFManager dof_manager;
    Index Nsublat = prim.basis.size();
    for(Index b = 0; b < Nsublat; b++) {
      if(prim.basis[b].site_occupant().size() > 1) {
        dof_manager.add_dof(prim.basis[b].site_occupant().type_name());
        break;
      }
    }
    for(Index b = 0; b < Nsublat; b++) {
      for(Index i = 0; i < prim.basis[i].displacement().size(); i++)
        dof_manager.add_dof(prim.basis[b].displacement()[i].type_name());
    }

    dof_manager.resize_neighborhood(nlist.size()*nlist.sublat_indices().size());

    // We can add more as needed
    dof_manager.register_dofs(tree);

    Index N_corr(tree.basis_set_size());
    std::stringstream private_def_stream, public_def_stream, interface_imp_stream, bfunc_imp_stream;

    std::string uclass_name;
    for(Index i = 0; i < class_name.size(); i++)
      uclass_name.push_back(std::toupper(class_name[i]));

    std::string indent(2, ' ');
    private_def_stream <<

                       indent << "  /// \\brief Clone the Clexulator\n" <<
                       indent << "  virtual " << class_name << "* _clone() const override {\n" <<
                       indent << "    return new " << class_name << "(*this);\n" <<
                       indent << "  }\n\n" <<

                       indent << "  // typedef for method pointers\n" <<
                       indent << "  typedef double (" << class_name << "::*BasisFuncPtr)() const;\n\n" <<

                       indent << "  // typedef for method pointers\n" <<
                       indent << "  typedef double (" << class_name << "::*DeltaBasisFuncPtr)(int, int) const;\n\n" <<

                       indent << "  // array of pointers to member functions for calculating basis functions\n" <<
                       indent << "  BasisFuncPtr m_orbit_func_list[" << N_corr << "];\n\n" <<

                       indent << "  // array of pointers to member functions for calculating flower functions\n" <<
                       indent << "  BasisFuncPtr m_flower_func_lists[" << Nsublat << "][" << N_corr << "];\n\n" <<

                       /**** for separate 1D method pointer lists:
                       indent << "  BasisFuncPtr
                       for(Index i = 0; i < Nsublat; i++) {
                         private_def_stream << " m_flower_func_at_" << i << "_list[" << N_corr << "]";
                         if(i + 1 < Nsublat)
                           private_def_stream << ',';
                       }
                       private_def_stream << ";\n\n" <<
                       **/

                       indent << "  // array of pointers to member functions for calculating DELTA flower functions\n" <<
                       indent << "  DeltaBasisFuncPtr m_delta_func_lists[" << Nsublat << "][" << N_corr << "];\n\n";

    /**** for separate 1D method pointer lists:
    indent << "  DeltaBasisFuncPtr";

    for(Index i = 0; i < Nsublat; i++) {
    private_def_stream << " m_delta_func_at_" << i << "_list[" << N_corr << "]";
    if(i + 1 < Nsublat)
    private_def_stream << ',';
    }
    private_def_stream << ";\n\n";
    **/

    dof_manager.print_clexulator_member_definitions(private_def_stream, tree, indent + "  ");

    // perhaps some more correlation calculating options here
    dof_manager.print_clexulator_private_method_definitions(private_def_stream, tree, indent + "  ");

    private_def_stream <<
                       indent << "  //default functions for basis function evaluation \n" <<
                       indent << "  double zero_func() const{ return 0.0;};\n" <<
                       indent << "  double zero_func(int,int) const{ return 0.0;};\n\n";

    public_def_stream <<
                      indent << "  " << class_name << "();\n\n" <<
                      indent << "  ~" << class_name << "();\n\n" <<

                      indent << "  /// \\brief Clone the " << class_name << "\n" <<
                      indent << "  std::unique_ptr<" << class_name << "> clone() const { \n" <<
                      indent << "    return std::unique_ptr<" << class_name << ">(_clone()); \n" <<
                      indent << "  }\n\n" <<

                      indent << "  /// \\brief Calculate contribution to global correlations from one unit cell\n" <<
                      indent << "  void calc_global_corr_contribution(double *corr_begin) const override;\n\n" <<

                      indent << "  /// \\brief Calculate contribution to select global correlations from one unit cell\n" <<
                      indent << "  void calc_restricted_global_corr_contribution(double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const override;\n\n" <<

                      indent << "  /// \\brief Calculate point correlations about basis site 'b_index'\n" <<
                      indent << "  void calc_point_corr(int b_index, double *corr_begin) const override;\n\n" <<

                      indent << "  /// \\brief Calculate select point correlations about basis site 'b_index'\n" <<
                      indent << "  void calc_restricted_point_corr(int b_index, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const override;\n\n" <<

                      indent << "  /// \\brief Calculate the change in point correlations due to changing an occupant\n" <<
                      indent << "  void calc_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin) const override;\n\n" <<

                      indent << "  /// \\brief Calculate the change in select point correlations due to changing an occupant\n" <<
                      indent << "  void calc_restricted_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const override;\n\n";

    dof_manager.print_clexulator_public_method_definitions(public_def_stream, tree, indent + "  ");


    //linear function index
    Index lf = 0, tlf;

    Array<FunctionVisitor *> labelers(dof_manager.get_function_label_visitors());
    //std::cout << "Initialized " << labelers.size() << " labelers \n";

    Array<std::string> orbit_method_names(N_corr);
    Array<Array<std::string> > flower_method_names(Nsublat, Array<std::string>(N_corr));
    //Array< Array<Array<std::string> > > dflower_method_names(N_corr, Array<Array<std::string> >(Nsublat));

    //this is very configuration-centric
    Array<Array<std::string> > dflower_method_names(Nsublat, Array<std::string>(N_corr));

    // temporary storage for formula
    Array<std::string> formulae, tformulae;

    bool make_newline(false);

    //loop over orbits
    for(Index np = 0; np < tree.size(); np++) {
      for(Index no = 0; no < tree[np].size(); no++) {
        if(np == 0)
          bfunc_imp_stream <<
                           indent << "// Basis functions for empty cluster:\n";
        else {
          bfunc_imp_stream <<
                           indent << "/**** Basis functions for orbit " << np << ", " << no << "****\n";
          tree[np][no].prototype.print(bfunc_imp_stream, '\n');
          bfunc_imp_stream << "****/\n";
        }

        formulae = tree[np][no].orbit_function_cpp_strings(labelers);
        tlf = formulae.size();

        make_newline = false;
        for(Index nf = 0; nf < formulae.size(); nf++) {
          if(!formulae[nf].size())
            continue;
          make_newline = true;
          orbit_method_names[lf + nf] = "eval_bfunc_" + std::to_string(np) + "_" + std::to_string(no) + "_" + std::to_string(nf);
          private_def_stream <<
                             indent << "  double " << orbit_method_names[lf + nf] << "() const;\n";

          bfunc_imp_stream <<
                           indent << "double " << class_name << "::" << orbit_method_names[lf + nf] << "() const{\n" <<
                           indent << "  return " << formulae[nf] << ";\n" <<
                           indent << "}\n";
        }
        if(make_newline) {
          bfunc_imp_stream << '\n';
          private_def_stream << '\n';
        }
        make_newline = false;

        // loop over flowers (i.e., basis sites of prim)
        const SiteOrbitBranch &asym_unit(tree.asym_unit());
        for(Index na = 0; na < asym_unit.size(); na++) {
          for(Index ne = 0; ne < asym_unit[na].size(); ne++) {
            Index nb = asym_unit[na][ne][0].basis_ind();
            auto nlist_index = find_index(nlist.sublat_indices(), nb);
            if(nlist_index != nlist.sublat_indices().size()) {
              formulae = tree[np][no].flower_function_cpp_strings(labelers, nlist_index);
              for(Index nf = 0; nf < formulae.size(); nf++) {
                if(!formulae[nf].size())
                  continue;
                make_newline = true;
                flower_method_names[nb][lf + nf] = "site_eval_at_" + std::to_string(nb) + "_bfunc_" + std::to_string(np) + "_" + std::to_string(no) + "_" + std::to_string(nf);
                private_def_stream <<
                                   indent << "  double " << flower_method_names[nb][lf + nf] << "() const;\n";

                bfunc_imp_stream <<
                                 indent << "double " << class_name << "::" << flower_method_names[nb][lf + nf] << "() const{\n" <<
                                 indent << "  return " << formulae[nf] << ";\n" <<
                                 indent << "}\n";

              }
            }
            if(make_newline) {
              bfunc_imp_stream << '\n';
              private_def_stream << '\n';
            }
            make_newline = false;

            // Very configuration-centric -> Find a way to move this block to OccupationDoFEnvironment:
            formulae.resize(formulae.size(), std::string());
            // loop over site basis functions
            const BasisSet &site_basis(asym_unit[na][ne].clust_basis);
            for(Index nsbf = 0; nsbf < site_basis.size(); nsbf++) {
              std::string delta_prefix = "(m_occ_func_" + std::to_string(nb) + "_" + std::to_string(nsbf) + "[occ_f] - m_occ_func_" + std::to_string(nb) + "_" + std::to_string(nsbf) + "[occ_i])";

              if(nlist_index != nlist.sublat_indices().size()) {
                tformulae = tree[np][no].delta_occfunc_flower_function_cpp_strings(site_basis, labelers, nlist_index, nb, nsbf);
                for(Index nf = 0; nf < tformulae.size(); nf++) {
                  if(!tformulae[nf].size())
                    continue;

                  if(formulae[nf].size())
                    formulae[nf] += " + ";

                  formulae[nf] += delta_prefix;

                  if(tformulae[nf] == "1" || tformulae[nf] == "(1)")
                    continue;

                  formulae[nf] += "*";
                  formulae[nf] += tformulae[nf];
                  //formulae[nf] += ")";
                }
              }
            }
            for(Index nf = 0; nf < formulae.size(); nf++) {
              if(!formulae[nf].size())
                continue;
              make_newline = true;

              dflower_method_names[nb][lf + nf] = "delta_site_eval_at_" + std::to_string(nb) + "_bfunc_" + std::to_string(np) + "_" + std::to_string(no) + "_" + std::to_string(nf);
              private_def_stream <<
                                 indent << "  double " << dflower_method_names[nb][lf + nf] << "(int occ_i, int occ_f) const;\n";

              bfunc_imp_stream <<
                               indent << "double " << class_name << "::" << dflower_method_names[nb][lf + nf] << "(int occ_i, int occ_f) const{\n" <<
                               indent << "  return " << formulae[nf] << ";\n" <<
                               indent << "}\n";
            }
            if(make_newline) {
              bfunc_imp_stream << '\n';
              private_def_stream << '\n';
            }
            make_newline = false;
            // \End Configuration specific part
          }
        }//\End loop over flowers

        lf += tlf;
      }
    }//Finished writing method definitions and implementations for basis functions

    //clean up:
    for(Index nl = 0; nl < labelers.size(); nl++)
      delete labelers[nl];
    labelers.clear();


    // Write constructor
    interface_imp_stream <<
                         indent << class_name << "::" << class_name << "() :\n" <<
                         indent << "  Clexulator_impl::Base(" << nlist.size() << ", " << N_corr << ") {\n";

    dof_manager.print_to_clexulator_constructor(interface_imp_stream, tree, indent + "  ");

    for(Index nf = 0; nf < orbit_method_names.size(); nf++) {
      if(orbit_method_names[nf].size() == 0)
        interface_imp_stream <<
                             indent << "  m_orbit_func_list[" << nf << "] = &" << class_name << "::zero_func;\n";
      else
        interface_imp_stream <<
                             indent << "  m_orbit_func_list[" << nf << "] = &" << class_name << "::" << orbit_method_names[nf] << ";\n";
    }
    interface_imp_stream << "\n\n";

    for(Index nb = 0; nb < flower_method_names.size(); nb++) {
      for(Index nf = 0; nf < flower_method_names[nb].size(); nf++) {
        if(flower_method_names[nb][nf].size() == 0)
          interface_imp_stream <<
                               indent << "  m_flower_func_lists[" << nb << "][" << nf << "] = &" << class_name << "::zero_func;\n";
        else
          interface_imp_stream <<
                               indent << "  m_flower_func_lists[" << nb << "][" << nf << "] = &" << class_name << "::" << flower_method_names[nb][nf] << ";\n";
      }
      interface_imp_stream << "\n\n";
    }

    for(Index nb = 0; nb < dflower_method_names.size(); nb++) {
      for(Index nf = 0; nf < dflower_method_names[nb].size(); nf++) {
        if(dflower_method_names[nb][nf].size() == 0)
          interface_imp_stream <<
                               indent << "  m_delta_func_lists[" << nb << "][" << nf << "] = &" << class_name << "::zero_func;\n";
        else
          interface_imp_stream <<
                               indent << "  m_delta_func_lists[" << nb << "][" << nf << "] = &" << class_name << "::" << dflower_method_names[nb][nf] << ";\n";
      }
      interface_imp_stream << "\n\n";
    }

    // Write weight matrix used for the neighbor list
    PrimNeighborList::Matrix3Type W = nlist.weight_matrix();
    interface_imp_stream << indent << "  m_weight_matrix.row(0) << " << W(0, 0) << ", " << W(0, 1) << ", " << W(0, 2) << ";\n";
    interface_imp_stream << indent << "  m_weight_matrix.row(1) << " << W(1, 0) << ", " << W(1, 1) << ", " << W(1, 2) << ";\n";
    interface_imp_stream << indent << "  m_weight_matrix.row(2) << " << W(2, 0) << ", " << W(2, 1) << ", " << W(2, 2) << ";\n\n";

    // Write neighborhood of UnitCellCoord
    // expand the nlist to contain 'global_orbitree' (all that is needed for now)
    std::set<UnitCellCoord> nbors;
    neighborhood(std::inserter(nbors, nbors.begin()), tree, prim, TOL);

    /*
    for(auto it = nbors.begin(); it != nbors.end(); ++it) {
      interface_imp_stream << indent << "  m_neighborhood.insert(UnitCellCoord("
                           << it->sublat() << ", "
                           << it->unitcell(0) << ", "
                           << it->unitcell(1) << ", "
                           << it->unitcell(2) << "));\n";
    }
    */

    interface_imp_stream << indent << "  m_neighborhood = std::set<UnitCellCoord> {\n";
    auto it = nbors.begin();
    while(it != nbors.end()) {
      interface_imp_stream << indent << "    {UnitCellCoord("
                           << it->sublat() << ", "
                           << it->unitcell(0) << ", "
                           << it->unitcell(1) << ", "
                           << it->unitcell(2) << ")}";
      ++it;
      if(it != nbors.end()) {
        interface_imp_stream << ",";
      }
      interface_imp_stream << "\n";
    }
    interface_imp_stream << indent << "  };\n";
    interface_imp_stream << "\n\n";

    interface_imp_stream << indent <<  "  m_orbit_neighborhood.resize(corr_size());\n";
    Index lno = 0;
    for(Index nb = 0; nb < tree.size(); ++nb) {
      for(Index no = 0; no < tree[nb].size(); ++no) {
        std::set<UnitCellCoord> orbit_nbors;
        orbit_neighborhood(std::inserter(orbit_nbors, orbit_nbors.begin()), tree, prim, nb, no, TOL);

        Index proto_index = lno;
        /*
        for(auto it = orbit_nbors.begin(); it != orbit_nbors.end(); ++it) {
          interface_imp_stream << indent << "  m_orbit_neighborhood[" << lno << "].insert(UnitCellCoord("
                               << it->sublat() << ", "
                               << it->unitcell(0) << ", "
                               << it->unitcell(1) << ", "
                               << it->unitcell(2) << "));\n";
        }
        */
        interface_imp_stream << indent << "  m_orbit_neighborhood[" << lno << "] = std::set<UnitCellCoord> {\n";
        auto it = orbit_nbors.begin();
        while(it != orbit_nbors.end()) {
          interface_imp_stream << indent << "    {UnitCellCoord("
                               << it->sublat() << ", "
                               << it->unitcell(0) << ", "
                               << it->unitcell(1) << ", "
                               << it->unitcell(2) << ")}";
          ++it;
          if(it != orbit_nbors.end()) {
            interface_imp_stream << ",";
          }
          interface_imp_stream << "\n";
        }
        interface_imp_stream << indent << "  };\n";
        ++lno;
        for(Index nf = 1; nf < tree.prototype(nb, no).clust_basis.size(); ++nf) {
          interface_imp_stream << indent << "  m_orbit_neighborhood[" << lno << "] = m_orbit_neighborhood[" << proto_index << "];\n";
          ++lno;
        }
        interface_imp_stream << "\n";

      }
    }


    interface_imp_stream <<
                         indent << "}\n\n";

    // Write destructor

    interface_imp_stream <<
                         indent << class_name << "::~" << class_name << "(){\n" <<

                         indent << "  //nothing here for now\n" <<

                         indent << "}\n\n";

    // Write evaluation methods

    interface_imp_stream <<
                         indent << "/// \\brief Calculate contribution to global correlations from one unit cell\n" <<
                         indent << "void " << class_name << "::calc_global_corr_contribution(double *corr_begin) const {\n" <<
                         indent << "  for(size_type i=0; i<corr_size(); i++){\n" <<
                         indent << "    *(corr_begin+i) = (this->*m_orbit_func_list[i])();\n" <<
                         indent << "  }\n" <<
                         indent << "}\n\n" <<

                         indent << "/// \\brief Calculate contribution to select global correlations from one unit cell\n" <<
                         indent << "void " << class_name << "::calc_restricted_global_corr_contribution(double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const {\n" <<
                         indent << "  for(; ind_list_begin<ind_list_end; ind_list_begin++){\n" <<
                         indent << "    *(corr_begin+*ind_list_begin) = (this->*m_orbit_func_list[*ind_list_begin])();\n" <<
                         indent << "  }\n" <<
                         indent << "}\n\n" <<

                         indent << "/// \\brief Calculate point correlations about basis site 'b_index'\n" <<
                         indent << "void " << class_name << "::calc_point_corr(int b_index, double *corr_begin) const {\n" <<
                         indent << "  for(size_type i=0; i<corr_size(); i++){\n" <<
                         indent << "    *(corr_begin+i) = (this->*m_flower_func_lists[b_index][i])();\n" <<
                         indent << "  }\n" <<
                         indent << "}\n\n" <<

                         indent << "/// \\brief Calculate select point correlations about basis site 'b_index'\n" <<
                         indent << "void " << class_name << "::calc_restricted_point_corr(int b_index, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const {\n" <<
                         indent << "  for(; ind_list_begin<ind_list_end; ind_list_begin++){\n" <<
                         indent << "    *(corr_begin+*ind_list_begin) = (this->*m_flower_func_lists[b_index][*ind_list_begin])();\n" <<
                         indent << "  }\n" <<
                         indent << "}\n\n" <<

                         indent << "/// \\brief Calculate the change in point correlations due to changing an occupant\n" <<
                         indent << "void " << class_name << "::calc_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin) const {\n" <<
                         indent << "  for(size_type i=0; i<corr_size(); i++){\n" <<
                         indent << "    *(corr_begin+i) = (this->*m_delta_func_lists[b_index][i])(occ_i, occ_f);\n" <<
                         indent << "  }\n" <<
                         indent << "}\n\n" <<

                         indent << "/// \\brief Calculate the change in select point correlations due to changing an occupant\n" <<
                         indent << "void " << class_name << "::calc_restricted_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const {\n" <<
                         indent << "  for(; ind_list_begin<ind_list_end; ind_list_begin++){\n" <<
                         indent << "    *(corr_begin+*ind_list_begin) = (this->*m_delta_func_lists[b_index][*ind_list_begin])(occ_i, occ_f);\n" <<
                         indent << "  }\n" <<
                         indent << "}\n\n";


    // PUT EVERYTHING TOGETHER
    stream <<
           "#include <cstddef>\n" <<
           "#include \"casm/clex/Clexulator.hh\"\n" <<
           "\n\n\n" <<
           "/****** CLEXULATOR CLASS FOR PRIM ******" << std::endl;

    jsonParser json;
    write_prim(prim, json, FRAC);
    stream << json;

    stream <<
           "**/\n\n\n" <<

           "/// \\brief Returns a Clexulator_impl::Base* owning a " << class_name << "\n" <<
           "extern \"C\" CASM::Clexulator_impl::Base* make_" + class_name << "();\n\n" <<

           "namespace CASM {\n\n" <<


           indent << "class " << class_name << " : public Clexulator_impl::Base {\n\n" <<

           indent << "public:\n\n" <<
           public_def_stream.str() << "\n" <<

           indent << "private:\n\n" <<
           private_def_stream.str() << "\n" <<

           indent << "};\n\n" << // close class definition

           indent <<

           "//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n" <<

           interface_imp_stream.str() <<
           bfunc_imp_stream.str() <<
           "}\n\n\n" <<      // close namespace

           "extern \"C\" {\n" <<
           indent << "/// \\brief Returns a Clexulator_impl::Base* owning a " << class_name << "\n" <<
           indent << "CASM::Clexulator_impl::Base* make_" + class_name << "() {\n" <<
           indent << "  return new CASM::" + class_name + "();\n" <<
           indent << "}\n\n" <<
           "}\n" <<

           "\n";
    // EOF

    return;
  }


}
