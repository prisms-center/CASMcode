#include "casm/clex/PrimClex.hh"

#include "casm/external/boost.hh"

#include "casm/misc/algorithm.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/ECIContainer.hh"
#include "casm/system/RuntimeLibrary.hh"
#include "casm/casm_io/SafeOfstream.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/app/AppIO.hh"

namespace CASM {
  //*******************************************************************************************
  //                                **** Constructors ****
  //*******************************************************************************************
  /// Initial construction of a PrimClex, from a primitive Structure
  PrimClex::PrimClex(const Structure &_prim, Log &log) :
    m_prim(_prim) {

    _init(log);

    return;
  }


  //*******************************************************************************************
  /// Construct PrimClex from existing CASM project directory
  ///  - read PrimClex and directory structure to generate all its Supercells and Configurations, etc.
  PrimClex::PrimClex(const fs::path &_root, Log &log):
    m_dir(_root),
    m_settings(_root),
    m_prim(read_prim(m_dir.prim())) {

    log.construct("CASM Project");
    log << "from: " << dir().root_dir() << "\n";

    _init(log);

  }

  /// Initialization routines
  ///  - If !root.empty(), read all saved data to generate all Supercells and Configurations, etc.
  void PrimClex::_init(Log &log) {

    std::vector<std::string> struc_mol_name = m_prim.get_struc_molecule_name();

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


    bool any_print = false;

    // here add stuff to read directory structure...

    // read param composition
    auto comp_axes = m_dir.composition_axes(m_settings.calctype(), m_settings.ref());
    if(fs::is_regular_file(comp_axes)) {
      log << "read: " << comp_axes << "\n";

      CompositionAxes opt(comp_axes);

      if(opt.has_current_axes) {
        m_has_composition_axes = true;
        m_comp_converter = opt.curr;
      }
    }

    // read chemical reference
    auto chem_ref_path = m_dir.chemical_reference(m_settings.calctype(), m_settings.ref());
    if(fs::is_regular_file(chem_ref_path)) {
      log << "read: " << chem_ref_path << "\n";
      m_chem_ref = notstd::make_cloneable<ChemicalReference>(read_chemical_reference(chem_ref_path, m_prim, settings().lin_alg_tol()));
    }

    // read supercells
    if(fs::is_regular_file(dir().SCEL())) {

      log << "read: " << dir().SCEL() << "\n";
      fs::ifstream scel(dir().SCEL());
      read_supercells(scel);

    }

    // read config_list
    if(fs::is_regular_file(dir().config_list())) {

      log << "read: " << dir().config_list() << "\n";
      read_config_list();
    }

    log << std::endl;
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
      std::cout << "Error in PrimClex::supercell(std::string scellname) const." << std::endl;
      std::cout << "  supercell '" << scellname << "' not found." << std::endl;
      exit(1);
    }
    return m_supercell_list[index];
  };

  //*******************************************************************************************
  /// Access supercell by name
  Supercell &PrimClex::supercell(std::string scellname) {
    Index index;
    if(!contains_supercell(scellname, index)) {
      std::cout << "Error in PrimClex::supercell(std::string scellname)." << std::endl;
      std::cout << "  supercell '" << scellname << "' not found." << std::endl;
      exit(1);
    }
    return m_supercell_list[index];
  };

  //*******************************************************************************************
  /// access configuration by name (of the form "scellname/[NUMBER]", e.g., ("SCEL1_1_1_1_0_0_0/0")
  const Configuration &PrimClex::configuration(const std::string &configname) const {
    std::vector<std::string> splt_vec;
    boost::split(splt_vec, configname, boost::is_any_of("/"), boost::token_compress_on);
    Index config_ind;
    if(splt_vec.size() != 2) {
      std::cerr << "CRITICAL ERROR: In PrimClex::configuration(), cannot locate configuration named '" << configname << "'\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }

    try {
      config_ind = boost::lexical_cast<Index>(splt_vec[1]);
    }
    catch(boost::bad_lexical_cast &) {
      std::cerr << "CRITICAL ERROR: In PrimClex::configuration(), malformed input:" << configname << "\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }


    return supercell(splt_vec[0]).config(config_ind);
  }

  //*******************************************************************************************

  Configuration &PrimClex::configuration(const std::string &configname) {
    std::vector<std::string> splt_vec;
    boost::split(splt_vec, configname, boost::is_any_of("/"), boost::token_compress_on);

    Index config_ind;
    if(splt_vec.size() != 2) {
      std::cerr << "CRITICAL ERROR: In PrimClex::configuration(), cannot locate configuration " << configname << "\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }

    try {
      config_ind = boost::lexical_cast<Index>(splt_vec[1]);
    }
    catch(boost::bad_lexical_cast &) {
      std::cerr << "CRITICAL ERROR: In PrimClex::configuration(), malformed input:" << configname << "\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }

    return supercell(splt_vec[0]).config(config_ind);
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


  /**
   * Generates all unique supercells between volStart and
   * volEnd of PRIM volumes. Call this routine before you
   * enumerate configurations.
   *
   * The provided transformation matrix can be used to enumerate
   * over lattice vectors that are not the ones belonging to the
   * primitive lattice (e.g. enumerate supercells of another supercell,
   * or enumerate over a particular lattice plane)
   *
   * The number of dimensions (must be equal to 1, 2 or 3) specify
   * which directions the supercells should be enumerated in, resulting
   * in 1D, 2D or 3D supercells. The enumeration is relative to the provided
   * transformation matrix. For dimension n, the first n columns of the
   * transformation matrix are used to construct supercells, while
   * the remaining 3-n columns remain fixed.
   *
   * The new functionality of restricted supercell enumeration can
   * be easily bypassed by passing dims=3 and G=Eigen::Matrix3i::Identity()
   *
   * @param[in] volStart Minimum volume supercell, relative to det(G)
   * @param[in] volEnd Maximum volume supercell, relative to det(G)
   * @param[in] dims Number of dimensions to enumerate over (1D, 2D or 3D supercells)
   * @param[in] G Generating matrix. Restricts enumeration to resulting vectors of P*G, where P=primitive.
   *
   */

  void PrimClex::generate_supercells(int volStart, int volEnd, int dims, const Eigen::Matrix3i &G, bool verbose) {
    Array < Lattice > supercell_lattices;
    m_prim.lattice().generate_supercells(supercell_lattices, m_prim.factor_group(), volStart, volEnd, dims, G);
    for(Index i = 0; i < supercell_lattices.size(); i++) {
      Index list_size = m_supercell_list.size();
      Index index = add_canonical_supercell(supercell_lattices[i]);
      if(m_supercell_list.size() != list_size && verbose) {
        std::cout << "  Generated: " << m_supercell_list[index].name() << "\n";
      }
      else {
        std::cout << "  Generated: " << m_supercell_list[index].name() << " (already existed)\n";
      }
    }
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
    return add_canonical_supercell(niggli(superlat, m_prim.point_group(), settings().crystallography_tol()));

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

  //*******************************************************************************************
  bool PrimClex::has_global_clexulator() const {
    if(!m_global_clexulator.initialized()) {
      if(!fs::exists(dir().clexulator_src(settings().name(), settings().bset()))) {
        return false;
      }
    }
    return true;
  }

  //*******************************************************************************************
  Clexulator PrimClex::global_clexulator(Log &status_log) const {
    if(!m_global_clexulator.initialized()) {

      if(!fs::exists(dir().clexulator_src(settings().name(), settings().bset()))) {
        throw std::runtime_error(
          std::string("Error loading clexulator ") + settings().bset() + ". No basis functions exist.");
      }

      m_global_clexulator = Clexulator(settings().global_clexulator(),
                                       dir().clexulator_dir(settings().bset()),
                                       nlist(),
                                       status_log,
                                       settings().compile_options(),
                                       settings().so_options());
    }
    return m_global_clexulator;
  }

  //*******************************************************************************************
  bool PrimClex::has_global_eci(std::string clex_name) const {

    if(m_global_eci.value().size()) {
      return true;
    }

    return fs::exists(dir().eci(clex_name,
                                settings().calctype(),
                                settings().ref(),
                                settings().bset(),
                                settings().eci()));
  }

  //*******************************************************************************************
  const ECIContainer &PrimClex::global_eci(std::string clex_name) const {
    if(!m_global_eci.value().size()) {
      fs::path eci_path = dir().eci(clex_name, settings().calctype(),
                                    settings().ref(), settings().bset(), settings().eci());
      if(!fs::exists(eci_path)) {
        throw std::runtime_error(
          std::string("Error loading global ECI. eci.json does not exist.\n")
          + "  Expected at: " + eci_path.string());
      }
      m_global_eci = read_eci(eci_path);
    }
    return m_global_eci;
  }

}

