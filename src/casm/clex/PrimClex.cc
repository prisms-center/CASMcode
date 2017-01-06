#include "casm/clex/PrimClex.hh"

#include "casm/external/boost.hh"

#include "casm/misc/algorithm.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/ECIContainer.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/system/RuntimeLibrary.hh"
#include "casm/misc/algorithm.hh"
#include "casm/casm_io/SafeOfstream.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/ECIContainer.hh"
#include "casm/app/AppIO.hh"

namespace CASM {
  //*******************************************************************************************
  //                                **** Constructors ****
  //*******************************************************************************************
  /// Initial construction of a PrimClex, from a primitive Structure
  PrimClex::PrimClex(const Structure &_prim, const Logging &logging) :
    Logging(logging),
    m_prim(_prim) {

    _init();

    return;
  }


  //*******************************************************************************************
  /// Construct PrimClex from existing CASM project directory
  ///  - read PrimClex and directory structure to generate all its Supercells and Configurations, etc.
  PrimClex::PrimClex(const fs::path &_root, const Logging &logging):
    Logging(logging),
    m_dir(_root),
    m_settings(_root),
    m_prim(read_prim(m_dir.prim())) {

    _init();

  }

  /// Initialization routines
  ///  - If !root.empty(), read all saved data to generate all Supercells and Configurations, etc.
  void PrimClex::_init() {

    log().construct("CASM Project");
    log() << "from: " << dir().root_dir() << "\n" << std::endl;

    auto struc_mol_name = prim().struc_molecule_name();
    m_vacancy_allowed = false;
    for(int i = 0; i < struc_mol_name.size(); ++i) {
      if(is_vacancy(struc_mol_name[i])) {
        m_vacancy_allowed = true;
        m_vacancy_index = i;
      }
    }

    if(dir().root_dir().empty()) {
      return;
    }

    bool read_settings = false;
    bool read_composition = true;
    bool read_chem_ref = true;
    bool read_configs = true;

    refresh(false, true, true, true);
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
        m_settings = ProjectSettings(dir().root_dir(), *this);
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
          m_chem_ref = notstd::make_cloneable<ChemicalReference>(read_chemical_reference(chem_ref_path, prim(), settings().lin_alg_tol()));
        }
      }
      catch(std::exception &e) {
        err_log().error("reading chemical_reference.json");
        err_log() << "file: " << chem_ref_path << "\n" << std::endl;
      }
    }

    if(read_configs) {

      m_supercell_list.clear();

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
        if(fs::is_regular_file(dir().config_list())) {
          log() << "read: " << dir().config_list() << "\n";
          read_config_list();
        }
      }
      catch(std::exception &e) {
        err_log().error("reading config_list.json");
        err_log() << "file: " << dir().config_list() << "\n" << std::endl;
      }
    }

    if(clear_clex) {
      m_nlist.reset();
      m_clex_basis.clear();
      m_clexulator.clear();
      m_eci.clear();
      log() << "refresh cluster expansions\n";
    }

    log() << std::endl;
  }


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
  const Structure &PrimClex::prim() const {
    return m_prim;
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
  boost::container::stable_vector<Supercell> &PrimClex::supercell_list() {
    return m_supercell_list;
  };

  //*******************************************************************************************
  /// const Access entire supercell_list
  const boost::container::stable_vector<Supercell> &PrimClex::supercell_list() const {
    return m_supercell_list;
  };

  //*******************************************************************************************
  /// const Access supercell by index
  const Supercell &PrimClex::supercell(Index i) const {
    return m_supercell_list[i];
  };

  //*******************************************************************************************
  /// Access supercell by index
  Supercell &PrimClex::supercell(Index i) {
    return m_supercell_list[i];
  };

  //*******************************************************************************************
  /// const Access supercell by name
  const Supercell &PrimClex::supercell(std::string scellname) const {
    Index index;
    if(!contains_supercell(scellname, index)) {
      err_log().error("Accessing supercell");
      err_log() << "supercell '" << scellname << "' not found." << std::endl;
      throw std::invalid_argument("Error in PrimClex::supercell(std::string scellname) const: Not found");
    }
    return m_supercell_list[index];
  };

  //*******************************************************************************************
  /// Access supercell by name
  Supercell &PrimClex::supercell(std::string scellname) {
    return const_cast<Supercell &>(static_cast<const PrimClex &>(*this).supercell(scellname));
  }

  //*******************************************************************************************
  /// Access supercell by Lattice, adding if necessary
  Supercell &PrimClex::supercell(const Lattice &lat) {
    return supercell(add_supercell(lat));
  }

  //*******************************************************************************************
  /// access configuration by name (of the form "scellname/[NUMBER]", e.g., ("SCEL1_1_1_1_0_0_0/0")
  const Configuration &PrimClex::configuration(const std::string &configname) const {
    auto res = Configuration::split_name(configname);

    try {
      return supercell(res.first).config_list().at(res.second);
    }
    catch(std::out_of_range &e) {
      err_log().error("Invalid config index");
      err_log() << "ERROR: In PrimClex::configuration(), configuration index out of range\n";
      err_log() << "configname: " << configname << "\n";
      err_log() << "index: " << res.second << "\n";
      err_log() << "config_list.size(): " << supercell(res.first).config_list().size() << "\n";
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
    if(m_supercell_list.size() == 0 || m_supercell_list[0].config_list().size() > 0)
      return config_iterator(this, 0, 0);
    return ++config_iterator(this, 0, 0);
  }

  //*******************************************************************************************
  /// Configuration iterator: end
  PrimClex::config_iterator PrimClex::config_end() {
    return config_iterator(this, m_supercell_list.size(), 0);
  }

  //*******************************************************************************************
  /// const Configuration iterator: begin
  PrimClex::config_const_iterator PrimClex::config_begin() const {
    if(m_supercell_list.size() == 0 || m_supercell_list[0].config_list().size() > 0)
      return config_const_iterator(this, 0, 0);
    return ++config_const_iterator(this, 0, 0);
  }

  //*******************************************************************************************
  /// const Configuration iterator: end
  PrimClex::config_const_iterator PrimClex::config_end() const {
    return config_const_iterator(this, m_supercell_list.size(), 0);
  }

  //*******************************************************************************************
  /// const Configuration iterator: begin
  PrimClex::config_const_iterator PrimClex::config_cbegin() const {
    if(m_supercell_list.size() == 0 || m_supercell_list[0].config_list().size() > 0)
      return config_const_iterator(this, 0, 0);
    return ++config_const_iterator(this, 0, 0);
  }

  //*******************************************************************************************
  /// const Configuration iterator: end
  PrimClex::config_const_iterator PrimClex::config_cend() const {
    return config_const_iterator(this, m_supercell_list.size(), 0);
  }

  //*******************************************************************************************
  /// Configuration iterator: begin
  PrimClex::config_iterator PrimClex::selected_config_begin() {
    //std::cout << "BEGINNING SELECTED CONFIG ITERATOR\n"
    //          << "m_supercell_list.size() is " << m_supercell_list.size() << "\n";

    if(m_supercell_list.size() == 0 || (m_supercell_list[0].config_list().size() > 0 && m_supercell_list[0].config(0).selected()))
      return config_iterator(this, 0, 0, true);
    return ++config_iterator(this, 0, 0, true);
  }

  //*******************************************************************************************
  /// Configuration iterator: end
  PrimClex::config_iterator PrimClex::selected_config_end() {
    return config_iterator(this, m_supercell_list.size(), 0, true);
  }

  //*******************************************************************************************
  /// const Configuration iterator: begin
  PrimClex::config_const_iterator PrimClex::selected_config_cbegin() const {
    if(m_supercell_list.size() == 0 || (m_supercell_list[0].config_list().size() > 0 && m_supercell_list[0].config(0).selected()))
      return config_const_iterator(this, 0, 0, true);
    return ++config_const_iterator(this, 0, 0, true);
  }

  //*******************************************************************************************
  /// const Configuration iterator: end
  PrimClex::config_const_iterator PrimClex::selected_config_cend() const {
    return config_const_iterator(this, m_supercell_list.size(), 0, true);
  }


  //*******************************************************************************************
  // **** IO ****
  //*******************************************************************************************
  /**
   * Re-write config_list.json, updating all the data
   */

  void PrimClex::write_config_list() {

    fs::path config_list_path = dir().config_list();
    if(m_supercell_list.size() == 0) {
      fs::remove(config_list_path);
      return;
    }

    jsonParser json;

    if(fs::exists(config_list_path)) {
      json.read(config_list_path);
    }
    else {
      json.put_obj();
    }

    for(Index s = 0; s < m_supercell_list.size(); s++) {
      m_supercell_list[s].write_config_list(json);
    }

    SafeOfstream file;
    file.open(config_list_path);
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
   *    return the index of the supercell in m_supercell_list
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
    if(!superlat.is_supercell_of(m_prim.lattice())) {
      std::cerr << "ERROR in PrimClex::add_canonical_supercell()." << std::endl
                << "  Passed a Supercell Lattice that is not a superlattice of PRIM lattice\n" << std::endl;
      assert(0);
      exit(1);
    }

    // temporary version, just check transf_mat equivalence
    // Does this check for equivalent supercells with different transformation matrices? Seems like it should
    // Insert second loop that goes over a symmetry operation list and applies it to the transformation matrix
    Supercell scel(this, superlat);
    scel.set_id(m_supercell_list.size());
    for(Index i = 0; i < m_supercell_list.size(); i++) {
      if(m_supercell_list[i].transf_mat() == scel.transf_mat())
        return i;
    }

    // if not already existing, add it
    m_supercell_list.push_back(scel);
    return m_supercell_list.size() - 1;
  }
  //*******************************************************************************************
  /**
   *  add a supercell if it doesn't already exist
   *  return the index of the supercell in m_supercell_list
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
    return add_canonical_supercell(canonical_equivalent_lattice(superlat, m_prim.point_group(), settings().crystallography_tol()));

  }

  //*******************************************************************************************
  void PrimClex::print_supercells() const {


    // also write supercells/supercell_path directories with LAT files
    try {
      fs::create_directory(dir().training_data());
    }
    catch(const fs::filesystem_error &ex) {
      std::cerr << "Error in PrimClex::print_supercells()." << std::endl;
      std::cerr << ex.what() << std::endl;
    }

    fs::ofstream scelfile(dir().SCEL());

    print_supercells(scelfile);

    for(Index i = 0; i < m_supercell_list.size(); i++) {
      try {
        fs::create_directory(dir().supercell_dir(m_supercell_list[i].name()));

        fs::path latpath = dir().LAT(m_supercell_list[i].name());
        if(!fs::exists(latpath)) {
          fs::ofstream latfile(latpath);
          m_supercell_list[i].real_super_lattice().print(latfile);
        }
      }
      catch(const fs::filesystem_error &ex) {
        std::cerr << "Error in PrimClex::print_supercells()." << std::endl;
        std::cerr << ex.what() << std::endl;
      }
    }
  }

  //*******************************************************************************************
  void PrimClex::print_supercells(std::ostream &stream) const {
    for(Index i = 0; i < m_supercell_list.size(); i++) {
      stream << "Supercell Name: '" << m_supercell_list[i].name() << "' Number: " << i << " Volume: " << m_supercell_list[i].transf_mat().determinant() << "\n";
      stream << "Supercell Transformation Matrix: \n";
      stream << m_supercell_list[i].transf_mat();
      stream << "\n";
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

    std::string s;
    while(!stream.eof()) {
      std::getline(stream, s);
      if(s[0] == 'S') {
        std::getline(stream, s);
        stream >> mat;

        add_canonical_supercell(Lattice(m_prim.lattice().lat_column_mat()*mat));
      }
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

    jsonParser json(dir().config_list());

    if(!json.contains("supercells")) {
      return;
    }

    for(Index i = 0; i < m_supercell_list.size(); i++) {
      m_supercell_list[i].read_config_list(json);
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
    for(Index s = 0; s < m_supercell_list.size(); s++) {
      amount_selected += m_supercell_list[s].amount_selected();
    }
    return amount_selected;
  }

  //*******************************************************************************************
  bool PrimClex::contains_supercell(std::string scellname, Index &index) const {
    for(Index i = 0; i < m_supercell_list.size(); i++) {
      if(m_supercell_list[i].name() == scellname) {
        index = i;
        return true;
      }
    }
    index = m_supercell_list.size();
    return false;

  };

  bool PrimClex::contains_supercell(const Supercell &scel) const {
    Index tmp;
    return contains_supercell(scel, tmp);
  }

  bool PrimClex::contains_supercell(const Supercell &scel, Index &index) const {
    return contains_supercell(scel.name(), index);
  }

  //*******************************************************************************************
  bool PrimClex::has_orbits(const ClexDescription &key) const {
    if(!fs::exists(dir().clust(key.bset))) {
      return false;
    }
    return true;
  }

  //*******************************************************************************************
  /// const Access to global orbitree
  bool PrimClex::has_clex_basis(const ClexDescription &key) const {
    auto it = m_clex_basis.find(key);
    if(it == m_clex_basis.end()) {
      if(!fs::exists(dir().clust(key.bset))) {
        return false;
      }
    }
    return true;

  };

  //*******************************************************************************************

  /// \brief Get iterators over the range of orbits
  const ClexBasis &PrimClex::clex_basis(const ClexDescription &key) const {

    auto it = m_clex_basis.find(key);
    if(it == m_clex_basis.end()) {

      it = m_clex_basis.insert(std::make_pair(key, ClexBasis(prim()))).first;

      std::vector<PrimPeriodicIntegralClusterOrbit> orbits;

      read_clust(
        std::back_inserter(orbits),
        jsonParser(dir().clust(key.bset)),
        prim(),
        prim().factor_group(),
        PrimPeriodicIntegralClusterSymCompare(settings().crystallography_tol()),
        settings().crystallography_tol()
      );

      jsonParser bspecs_json;
      bspecs_json.read(dir().bspecs(key.bset));

      ClexBasis &clex_basis = it->second;
      clex_basis.generate(orbits.begin(), orbits.end(), bspecs_json);

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

      it = m_clexulator.insert(
             std::make_pair(key, Clexulator(settings().name() + "_Clexulator",
                                            dir().clexulator_dir(key.bset),
                                            nlist(),
                                            log(),
                                            settings().compile_options(),
                                            settings().so_options()))).first;
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

}

