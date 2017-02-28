#include "casm/clex/PrimClex.hh"

#include "casm/external/boost.hh"

#include "casm/clex/ECIContainer.hh"
#include "casm/system/RuntimeLibrary.hh"
#include "casm/casm_io/SafeOfstream.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clex/ECIContainer.hh"

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

      if(!m_db_handler) {
        m_db_handler = notstd::make_unique<DB::DatabaseHandler>(*this);
      }
      else {
        // lazy initialization means we just need to close, and the db will be
        // re-opened when needed
        m_db_handler->close();
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
      clex_basis.generate(orbits.begin(), orbits.end(), bspecs_json, {"occupation"});

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

