#include "casm/clex/PrimClex_impl.hh"

#include "casm/casm_io/SafeOfstream.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/ChemicalReference.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/ECIContainer.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/database/DatabaseHandler_impl.hh"
#include "casm/database/DatabaseTypes_impl.hh"

#include "casm/app/AppIO.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ClexDescription.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/EnumeratorHandler_impl.hh"
#include "casm/app/QueryHandler_impl.hh"

namespace CASM {

  struct PrimClex::PrimClexData {

    typedef PrimClex::PrimType PrimType;

    PrimClexData(const Structure &_prim) :
      prim_ptr(std::make_shared<PrimType>(_prim)) {
      //Guarantee presence of symmetry info;
      prim_ptr->factor_group();
    }

    PrimClexData(const fs::path &_root) :
      dir(_root),
      settings(_root),
      prim_ptr(std::make_shared<PrimType>(read_prim(dir.prim(), settings.crystallography_tol(), &(settings.hamiltonian_modules())))) {

      //Guarantee presence of symmetry info;
      prim_ptr->factor_group();

    }

    ~PrimClexData() {}

    DirectoryStructure dir;
    ProjectSettings settings;

    std::shared_ptr<const PrimType> prim_ptr;
    bool vacancy_allowed;
    Index vacancy_index;

    std::unique_ptr<DB::DatabaseHandler> db_handler;

    /// CompositionConverter specifies parameteric composition axes and converts between
    ///   parametric composition and mol composition
    bool has_composition_axes = false;
    CompositionConverter comp_converter;

    /// ChemicalReference specifies a reference for formation energies, chemical
    /// potentials, etc.
    notstd::cloneable_ptr<ChemicalReference> chem_ref;

    /// Stores the neighboring UnitCell and which sublattices to include in neighbor lists
    /// - mutable for lazy construction
    mutable notstd::cloneable_ptr<PrimNeighborList> nlist;

    mutable std::map<ClexDescription, ClexBasis> clex_basis;
    mutable std::map<ClexDescription, Clexulator> clexulator;
    mutable std::map<ClexDescription, ECIContainer> eci;

  };

  //  **** Constructors ****

  /// Initial construction of a PrimClex, from a primitive Structure
  PrimClex::PrimClex(const Structure &_prim, const Logging &logging) :
    Logging(logging),
    m_data(new PrimClexData(_prim)) {

    m_data->settings.set_crystallography_tol(TOL);

    _init();

    return;
  }

  /// Construct PrimClex from existing CASM project directory
  ///  - read PrimClex and directory structure to generate all its Supercells and Configurations, etc.
  PrimClex::PrimClex(const fs::path &_root, const Logging &logging):
    Logging(logging),
    m_data(new PrimClexData(_root)) {

    _init();

  }

  /// Necessary for "pointer to implementation"
  PrimClex::~PrimClex() {}

  /// Initialization routines
  ///  - If !root.empty(), read all saved data to generate all Supercells and Configurations, etc.
  void PrimClex::_init() {

    log().construct("CASM Project");
    log() << "from: " << dir().root_dir() << "\n" << std::endl;

    auto struc_mol_name = struc_molecule_name(prim());
    m_data->vacancy_allowed = false;
    for(int i = 0; i < struc_mol_name.size(); ++i) {
      if(xtal::is_vacancy(struc_mol_name[i])) {
        m_data->vacancy_allowed = true;
        m_data->vacancy_index = i;
      }
    }

    if(dir().root_dir().empty()) {
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
        m_data->settings = ProjectSettings(dir().root_dir(), *this);
      }
      catch(std::exception &e) {
        err_log().error("reading project_settings.json");
        err_log() << "file: " << m_data->dir.project_settings() << "\n" << std::endl;
      }
    }

    if(read_composition) {
      m_data->has_composition_axes = false;
      auto comp_axes = m_data->dir.composition_axes();

      try {
        if(fs::is_regular_file(comp_axes)) {
          log() << "read: " << comp_axes << "\n";

          CompositionAxes opt(comp_axes);

          if(opt.has_current_axes) {
            m_data->has_composition_axes = true;
            m_data->comp_converter = opt.curr;
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
      m_data->chem_ref.reset();
      auto chem_ref_path = m_data->dir.chemical_reference(m_data->settings.default_clex().calctype, m_data->settings.default_clex().ref);

      try {
        if(fs::is_regular_file(chem_ref_path)) {
          log() << "read: " << chem_ref_path << "\n";
          m_data->chem_ref = notstd::make_cloneable<ChemicalReference>(read_chemical_reference(chem_ref_path, prim(), settings().lin_alg_tol()));
        }
      }
      catch(std::exception &e) {
        err_log().error("reading chemical_reference.json");
        err_log() << "file: " << chem_ref_path << "\n" << std::endl;
      }
    }

    if(read_configs) {

      if(!m_data->db_handler) {
        m_data->db_handler = notstd::make_unique<DB::DatabaseHandler>(*this);
      }
      else {
        // lazy initialization means we just need to close, and the db will be
        // re-opened when needed
        m_data->db_handler->close();
      }
    }

    if(clear_clex) {
      m_data->nlist.reset();
      m_data->clex_basis.clear();
      m_data->clexulator.clear();
      m_data->eci.clear();
      log() << "refresh cluster expansions\n";
    }

    log() << std::endl;
  }

  const DirectoryStructure &PrimClex::dir() const {
    return m_data->dir;
  }

  ProjectSettings &PrimClex::settings() {
    return m_data->settings;
  }

  const ProjectSettings &PrimClex::settings() const {
    return m_data->settings;
  }

  /// \brief Get the crystallography_tol
  double PrimClex::crystallography_tol() const {
    return prim().lattice().tol();
  }

  // ** Composition accessors **

  /// const Access CompositionConverter object
  bool PrimClex::has_composition_axes() const {
    return m_data->has_composition_axes;
  }

  /// const Access CompositionConverter object
  const CompositionConverter &PrimClex::composition_axes() const {
    return m_data->comp_converter;
  }

  // ** Chemical reference **

  /// check if ChemicalReference object initialized
  bool PrimClex::has_chemical_reference() const {
    return static_cast<bool>(m_data->chem_ref);
  }

  /// const Access ChemicalReference object
  const ChemicalReference &PrimClex::chemical_reference() const {
    return *m_data->chem_ref;
  }


  // ** Prim and Orbitree accessors **

  /// const Access to primitive Structure
  const PrimClex::PrimType &PrimClex::prim() const {
    return *(m_data->prim_ptr);
  }

  std::shared_ptr<const PrimClex::PrimType> &PrimClex::shared_prim() const {
    return this->m_data->prim_ptr;
  }

  Index PrimClex::n_basis() const {
    return prim().basis().size();
  }

  PrimNeighborList &PrimClex::nlist() const {

    // lazy neighbor list generation
    if(!m_data->nlist) {

      // construct nlist
      m_data->nlist = notstd::make_cloneable<PrimNeighborList>(
                        settings().nlist_weight_matrix(),
                        settings().nlist_sublat_indices().begin(),
                        settings().nlist_sublat_indices().end()
                      );
    }

    return *m_data->nlist;
  }

  /// returns true if vacancy are an allowed species
  bool PrimClex::vacancy_allowed() const {
    return m_data->vacancy_allowed;
  }

  /// returns the index of vacancies in composition vectors
  Index PrimClex::vacancy_index() const {
    return m_data->vacancy_index;
  }

  template<typename T>
  DB::ValDatabase<T> &PrimClex::generic_db() const {
    return m_data->db_handler->template generic_db<T>();
  }

  template<typename T>
  const DB::ValDatabase<T> &PrimClex::const_generic_db() const {
    return m_data->db_handler->template const_generic_db<T>();
  }

  template<typename T>
  DB::Database<T> &PrimClex::db() const {
    return m_data->db_handler->template db<T>();
  }

  template<typename T>
  const DB::Database<T> &PrimClex::const_db() const {
    return m_data->db_handler->template const_db<T>();
  }

  template<typename T>
  DB::PropertiesDatabase &PrimClex::db_props(std::string calc_type) const {
    return m_data->db_handler->template db_props<T>(calc_type);
  }

  template<typename T>
  const DB::PropertiesDatabase &PrimClex::const_db_props(std::string calc_type) const {
    return m_data->db_handler->template const_db_props<T>(calc_type);
  }

  DB::DatabaseHandler &PrimClex::db_handler() const {
    return *m_data->db_handler;
  }

  const DB::DatabaseHandler &PrimClex::const_db_handler() const {
    return *m_data->db_handler;
  }

  bool PrimClex::has_orbits(const ClexDescription &key) const {
    if(!fs::exists(dir().clust(key.bset))) {
      return false;
    }
    return true;
  }

  /// const Access to global orbitree
  bool PrimClex::has_clex_basis(const ClexDescription &key) const {
    auto it = m_data->clex_basis.find(key);
    if(it == m_data->clex_basis.end()) {
      if(!fs::exists(dir().clust(key.bset))) {
        return false;
      }
    }
    return true;

  };

  /// \brief Get iterators over the range of orbits
  const ClexBasis &PrimClex::clex_basis(const ClexDescription &key) const {

    auto it = m_data->clex_basis.find(key);
    if(it == m_data->clex_basis.end()) {

      jsonParser bspecs_json;
      bspecs_json.read(dir().bspecs(key.bset));

      it = m_data->clex_basis.insert(std::make_pair(key, ClexBasis(prim(), bspecs_json))).first;

      std::vector<PrimPeriodicOrbit<IntegralCluster>> orbits;

      typedef PrimPeriodicSymCompare<IntegralCluster> symcompare_type;

      read_clust(
        std::back_inserter(orbits),
        jsonParser(dir().clust(key.bset)),
        prim(),
        prim().factor_group(),
        symcompare_type(this->shared_prim(), crystallography_tol()),
        crystallography_tol()
      );

      ClexBasis &clex_basis = it->second;
      clex_basis.generate(orbits.begin(), orbits.end(), bspecs_json);

    }

    return it->second;

  }

  bool PrimClex::has_clexulator(const ClexDescription &key) const {
    auto it = m_data->clexulator.find(key);
    if(it == m_data->clexulator.end()) {
      if(!fs::exists(dir().clexulator_src(settings().name(), key.bset))) {
        return false;
      }
    }
    return true;
  }

  Clexulator PrimClex::clexulator(const ClexDescription &key) const {

    auto it = m_data->clexulator.find(key);
    if(it == m_data->clexulator.end()) {

      if(!fs::exists(dir().clexulator_src(settings().name(), key.bset))) {
        throw std::runtime_error(
          std::string("Error loading clexulator ") + key.bset + ". No basis functions exist.");
      }

      try {
        it = m_data->clexulator.insert(
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

  bool PrimClex::has_eci(const ClexDescription &key) const {

    auto it = m_data->eci.find(key);
    if(it == m_data->eci.end()) {
      return fs::exists(dir().eci(key.property, key.calctype, key.ref, key.bset, key.eci));
    }
    return true;
  }

  const ECIContainer &PrimClex::eci(const ClexDescription &key) const {

    auto it = m_data->eci.find(key);
    if(it == m_data->eci.end()) {
      fs::path eci_path = dir().eci(key.property, key.calctype, key.ref, key.bset, key.eci);
      if(!fs::exists(eci_path)) {
        throw std::runtime_error(
          std::string("Error loading ECI. eci.json does not exist.\n")
          + "  Expected at: " + eci_path.string());
      }

      it = m_data->eci.insert(std::make_pair(key, read_eci(eci_path))).first;
    }
    return it->second;
  }
}

#define INST_PrimClex_orbits_vec(OrbitType, SymCompareType) \
  template std::back_insert_iterator<std::vector<OrbitType>> PrimClex::orbits( \
    const ClexDescription &, \
    std::back_insert_iterator<std::vector<OrbitType>>, \
    const SymCompareType&) const;

namespace CASM {
  INST_PrimClex_orbits_vec(PrimPeriodicOrbit<IntegralCluster>, PrimPeriodicSymCompare<IntegralCluster>);
}

#include "casm/database/DatabaseTypes.hh"

// explicit template instantiations
#define INST_PrimClex(r, data, type) \
template DB::Database<type> &PrimClex::db<type>() const; \
template const DB::Database<type> &PrimClex::const_db<type>() const; \
template DB::ValDatabase<type> &PrimClex::generic_db<type>() const; \
template const DB::ValDatabase<type> &PrimClex::const_generic_db<type>() const; \

// explicit template instantiations
#define INST_PrimClexProps(r, data, type) \
template DB::PropertiesDatabase &PrimClex::db_props<type>(std::string calc_type) const; \
template const DB::PropertiesDatabase &PrimClex::const_db_props<type>(std::string calc_type) const; \

namespace CASM {
  BOOST_PP_SEQ_FOR_EACH(INST_PrimClex, _, CASM_DB_TYPES)
  BOOST_PP_SEQ_FOR_EACH(INST_PrimClexProps, _, CASM_DB_CONFIG_TYPES)
}
