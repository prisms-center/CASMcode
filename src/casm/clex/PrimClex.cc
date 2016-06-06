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

    std::vector<std::string> struc_mol_name = prim().get_struc_molecule_name();

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
      m_chem_ref = notstd::make_cloneable<ChemicalReference>(read_chemical_reference(chem_ref_path, prim, settings().lin_alg_tol()));
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
    prim.lattice().generate_supercells(supercell_lattices, prim.factor_group(), volStart, volEnd, dims, G);
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
    return add_canonical_supercell(niggli(superlat, prim.point_group(), crystallography_tol()));

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

        add_canonical_supercell(Lattice(prim.lattice().lat_column_mat()*mat));
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
  Clexulator PrimClex::global_clexulator() const {
    if(!m_global_clexulator.initialized()) {

      if(!fs::exists(dir().clexulator_src(settings().name(), settings().bset()))) {
        throw std::runtime_error(
          std::string("Error loading clexulator ") + settings().bset() + ". No basis functions exist.");
      }

      m_global_clexulator = Clexulator(settings().global_clexulator(),
                                       dir().clexulator_dir(settings().bset()),
                                       nlist(),
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

  //*******************************************************************************************
  /// \brief Make orbitree. For now specifically global.
  ///
  /// \param prim Primitive Structure. non-const due to Structure::set_site_internals.
  /// \param json bspecs.json file
  SiteOrbitree make_orbitree(Structure &prim, const jsonParser &json) {

    try {

      SiteOrbitree tree(prim.lattice());
      tree.set_bspecs(json);
      // --- first generate global orbitree -------------


      //global_orbitree.read_CSPECS(in_clust);
      tree.min_num_components = 2;
      tree.min_length = CASM::TOL;

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

    for(auto it = nbors.begin(); it != nbors.end(); ++it) {
      interface_imp_stream << indent << "  m_neighborhood.insert(UnitCellCoord("
                           << it->sublat() << ", "
                           << it->unitcell(0) << ", "
                           << it->unitcell(1) << ", "
                           << it->unitcell(2) << "));\n";
    }
    interface_imp_stream << "\n\n";

    interface_imp_stream << indent <<  "  m_orbit_neighborhood.resize(corr_size());\n";
    Index lno = 0;
    for(Index nb = 0; nb < tree.size(); ++nb) {
      for(Index no = 0; no < tree[nb].size(); ++no) {
        std::set<UnitCellCoord> orbit_nbors;
        orbit_neighborhood(std::inserter(orbit_nbors, orbit_nbors.begin()), tree, prim, nb, no, TOL);

        Index proto_index = lno;
        for(auto it = orbit_nbors.begin(); it != orbit_nbors.end(); ++it) {
          interface_imp_stream << indent << "  m_orbit_neighborhood[" << lno << "].insert(UnitCellCoord("
                               << it->sublat() << ", "
                               << it->unitcell(0) << ", "
                               << it->unitcell(1) << ", "
                               << it->unitcell(2) << "));\n";
        }
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

