#include "casm/clex/PrimClex_impl.hh"

#include "casm/casm_io/SafeOfstream.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/ChemicalReference.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/ClexBasisWriter_impl.hh"
#include "casm/clex/ECIContainer.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clusterography/IntegralCluster_impl.hh"
#include "casm/database/DatabaseHandler_impl.hh"
#include "casm/database/DatabaseTypes_impl.hh"

#include "casm/app/AppIO_impl.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ClexDescription.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/EnumeratorHandler_impl.hh"
#include "casm/app/HamiltonianModules_impl.hh"
#include "casm/app/QueryHandler_impl.hh"
#include <memory>

namespace CASM {

  BasicStructure read_prim(ProjectSettings const &project_settings) {
    return read_prim(
             project_settings.dir().prim(),
             project_settings.crystallography_tol(),
             &project_settings.hamiltonian_modules());
  }

  std::shared_ptr<Structure const> read_shared_prim(ProjectSettings const &project_settings) {
    return std::make_shared<Structure const>(read_prim(project_settings));
  }


  struct PrimClex::PrimClexData {

    typedef PrimClex::PrimType PrimType;
    typedef std::shared_ptr<PrimType const> PrimType_ptr;

    PrimClexData(
      ProjectSettings const &_project_settings,
      std::shared_ptr<PrimType const> _shared_prim) :
      settings(_project_settings),
      prim_ptr(_shared_prim) {
      //Guarantee presence of symmetry info;
      prim_ptr->factor_group();
    }

    PrimClexData(const fs::path &_root) :
      settings(open_project_settings(_root)),
      prim_ptr(read_shared_prim(settings)) {

      //Guarantee presence of symmetry info;
      prim_ptr->factor_group();

    }

    ~PrimClexData() {}

    ProjectSettings settings;

    PrimType_ptr prim_ptr;
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

    typedef std::string BasisSetName;
    mutable std::map<BasisSetName, jsonParser> basis_set_specs;
    mutable std::map<BasisSetName, ClexBasis> clex_basis;
    mutable std::map<BasisSetName, Clexulator> clexulator;
    mutable std::map<ClexDescription, ECIContainer> eci;

  };

  //  **** Constructors ****

  /// Initial construction of a PrimClex, from a primitive Structure
  PrimClex::PrimClex(
    ProjectSettings const &_project_settings,
    std::shared_ptr<PrimType const> _shared_prim,
    const Logging &logging) :
    Logging(logging),
    m_data(new PrimClexData(_project_settings, _shared_prim)) {

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

    if(has_dir()) {
      log().construct("CASM Project");
      log() << "from: " << dir().root_dir() << "\n" << std::endl;
    }

    auto struc_mol_name = xtal::struc_molecule_name(prim());
    m_data->vacancy_allowed = false;
    for(int i = 0; i < struc_mol_name.size(); ++i) {
      if(xtal::is_vacancy(struc_mol_name[i])) {
        m_data->vacancy_allowed = true;
        m_data->vacancy_index = i;
      }
    }

    if(!has_dir()) {
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
        m_data->settings = open_project_settings(dir().root_dir());
      }
      catch(std::exception &e) {
        err_log().error("reading project_settings.json");
        err_log() << "file: " << dir().project_settings() << "\n" << std::endl;
        throw e;
      }
    }

    if(read_composition) {
      m_data->has_composition_axes = false;
      auto comp_axes = dir().composition_axes();

      try {
        if(fs::is_regular_file(comp_axes)) {
          log() << "read: " << comp_axes << "\n";

          CompositionAxes opt(comp_axes);

          if(opt.has_current_axes()) {
            m_data->has_composition_axes = true;
            m_data->comp_converter = opt.curr;
          }
        }
      }
      catch(std::exception &e) {
        err_log().error("reading composition_axes.json");
        err_log() << "file: " << comp_axes << "\n" << std::endl;
        throw e;
      }
    }

    if(read_chem_ref) {

      // read chemical reference
      m_data->chem_ref.reset();
      auto chem_ref_path = dir().chemical_reference(m_data->settings.default_clex().calctype, m_data->settings.default_clex().ref);

      try {
        if(fs::is_regular_file(chem_ref_path)) {
          log() << "read: " << chem_ref_path << "\n";
          m_data->chem_ref = notstd::make_cloneable<ChemicalReference>(read_chemical_reference(chem_ref_path, prim(), settings().lin_alg_tol()));
        }
      }
      catch(std::exception &e) {
        err_log().error("reading chemical_reference.json");
        err_log() << "file: " << chem_ref_path << "\n" << std::endl;
        throw e;
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

  ProjectSettings &PrimClex::settings() {
    return m_data->settings;
  }

  const ProjectSettings &PrimClex::settings() const {
    return m_data->settings;
  }

  bool PrimClex::has_dir() const {
    return settings().has_dir();
  }

  const DirectoryStructure &PrimClex::dir() const {
    return settings().dir();
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

  std::shared_ptr<PrimClex::PrimType const> &PrimClex::shared_prim() const {
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

  bool PrimClex::has_basis_set_specs(std::string const &basis_set_name) const {
    auto it = m_data->basis_set_specs.find(basis_set_name);
    if(it == m_data->basis_set_specs.end()) {
      if(!fs::exists(dir().bspecs(basis_set_name))) {
        return false;
      }
    }
    return true;
  }

  jsonParser const &PrimClex::basis_set_specs(std::string const &basis_set_name) const {
    auto it = m_data->basis_set_specs.find(basis_set_name);
    if(it == m_data->basis_set_specs.end()) {
      if(!fs::exists(dir().bspecs(basis_set_name))) {
        std::stringstream ss;
        ss << "Error in PrimClex::basis_set_specs: basis set named '" << basis_set_name
           << "' does not exist.";
        throw std::runtime_error(ss.str());
      }

      fs::path basis_set_specs_path = dir().bspecs(basis_set_name);
      try {
        jsonParser basis_set_specs {basis_set_specs_path};
        it = m_data->basis_set_specs.emplace(basis_set_name, basis_set_specs).first;
      }
      catch(std::exception &e) {
        err_log().error("reading bspecs.json");
        err_log() << "file: " << basis_set_specs_path << "\n" << std::endl;
        throw e;
      }
    }
    return it->second;
  }

  bool PrimClex::has_orbits(std::string const &basis_set_name) const {
    if(!fs::exists(dir().clust(basis_set_name))) {
      return false;
    }
    return true;
  }

  /// const Access to global orbitree
  bool PrimClex::has_clex_basis(std::string const &basis_set_name) const {
    auto it = m_data->clex_basis.find(basis_set_name);
    if(it == m_data->clex_basis.end()) {
      if(!fs::exists(dir().clust(basis_set_name))) {
        return false;
      }
    }
    return true;

  };

  /// \brief Get iterators over the range of orbits
  const ClexBasis &PrimClex::clex_basis(std::string const &basis_set_name) const {

    auto it = m_data->clex_basis.find(basis_set_name);
    if(it == m_data->clex_basis.end()) {

      jsonParser basis_set_specs {dir().bspecs(basis_set_name)};
      it = m_data->clex_basis.emplace(
             basis_set_name, ClexBasis {this->shared_prim(), basis_set_specs}).first;

      std::vector<PrimPeriodicOrbit<IntegralCluster>> orbits;

      typedef PrimPeriodicSymCompare<IntegralCluster> symcompare_type;

      read_clust(
        std::back_inserter(orbits),
        jsonParser(dir().clust(basis_set_name)),
        prim(),
        prim().factor_group(),
        symcompare_type(this->shared_prim(), crystallography_tol()),
        crystallography_tol()
      );

      ClexBasis &clex_basis = it->second;
      clex_basis.generate(orbits.begin(), orbits.end(), basis_set_specs);

    }

    return it->second;

  }

  bool PrimClex::has_clexulator(std::string const &basis_set_name) const {
    auto it = m_data->clexulator.find(basis_set_name);
    if(it == m_data->clexulator.end()) {
      if(!fs::exists(dir().clexulator_src(settings().project_name(), basis_set_name))) {
        return false;
      }
    }
    return true;
  }

  Clexulator PrimClex::clexulator(std::string const &basis_set_name) const {

    auto it = m_data->clexulator.find(basis_set_name);
    if(it == m_data->clexulator.end()) {

      if(!fs::exists(dir().clexulator_src(settings().project_name(), basis_set_name))) {
        throw std::runtime_error(
          std::string("Error loading clexulator ") + basis_set_name + ". No basis functions exist.");
      }

      try {
        Clexulator clexulator = read_clexulator(settings(), basis_set_name, nlist());
        it = m_data->clexulator.emplace(basis_set_name, clexulator).first;
      }
      catch(std::exception &e) {
        // TODO: not sure why this fails...
        // log() << "Error constructing Clexulator. Current settings: \n" << std::endl;
        // settings().print_compiler_settings_summary(log());

        std::cout << "Error constructing Clexulator. Current settings: \n" << std::endl;
        Log tlog(std::cout);
        print_compiler_settings_summary(settings(), tlog);
        throw e;
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

  void _throw_if_no_bset(std::string const &basis_set_name, DirectoryStructure const &dir) {
    if(dir.root_dir().empty()) {
      throw std::runtime_error("Error accessing bset." + basis_set_name +
                               ": No root directory set for project.");
    }
    auto basis_set_specs_path = dir.bspecs(basis_set_name);
    if(!fs::exists(basis_set_specs_path)) {
      throw std::runtime_error("Error accessing bset." + basis_set_name +
                               ": Does not exist." + "  Checked for bspecs at: " + basis_set_specs_path.string());
    }
  }

  template<typename OrbitVecType>
  void _write_clexulator(
    std::shared_ptr<Structure const> shared_prim,
    ProjectSettings const &settings,
    std::string const &basis_set_name,
    jsonParser const &basis_set_specs,
    PrimNeighborList &prim_neighbor_list,
    OrbitVecType const &orbits,
    ClexBasis const &clex_basis) {

    const auto &dir = settings.dir();
    double xtal_tol = settings.crystallography_tol();
    Log &log = CASM::log();

    // write clust
    jsonParser clust_json;
    write_clust(orbits.begin(), orbits.end(), clust_json, ProtoSitesPrinter(), basis_set_specs);
    clust_json.write(dir.clust(basis_set_name));
    log.write(dir.clust(basis_set_name).string());
    log << std::endl;

    // write basis
    jsonParser basis_json;
    write_site_basis_funcs(shared_prim, clex_basis, basis_json);
    write_clust(
      orbits.begin(),
      orbits.end(),
      basis_json,
      ProtoFuncsPrinter {clex_basis, shared_prim->shared_structure()},
      basis_set_specs);
    basis_json.write(dir.basis(basis_set_name));
    log.write(dir.basis(basis_set_name).string());
    log << std::endl;

    // write source code
    fs::ofstream outfile;
    outfile.open(dir.clexulator_src(settings.project_name(), basis_set_name));
    std::string parampack_type {"DEFAULT"};
    basis_set_specs.get_if(parampack_type, "param_pack");
    ClexBasisWriter clexwriter {*shared_prim, parampack_type};
    clexwriter.print_clexulator(
      settings.global_clexulator_name(), clex_basis, orbits, prim_neighbor_list, outfile, xtal_tol);
    outfile.close();
    log.write(dir.clexulator_src(settings.project_name(), basis_set_name).string());
    log << std::endl;
  }

  void write_clexulator(
    std::shared_ptr<Structure const> shared_prim,
    ProjectSettings const &settings,
    std::string const &basis_set_name,
    jsonParser const &basis_set_specs,
    PrimNeighborList &prim_neighbor_list) {

    auto const &dir = settings.dir();
    _throw_if_no_bset(basis_set_name, dir);
    dir.delete_bset_data(settings.project_name(), basis_set_name);

    // TODO: update this function with ClusterSpecs

    CLUSTER_PERIODICITY_TYPE clex_periodicity_type = basis_set_specs.contains("local_bpsecs") ?
                                                     CLUSTER_PERIODICITY_TYPE::LOCAL : CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC;

    if(clex_periodicity_type == CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC) {
      std::vector<PrimPeriodicIntegralClusterOrbit> orbits;

      log().construct("Orbitree");
      log() << std::endl;

      std::vector<DoFKey> _dofs;
      if(basis_set_specs.contains("basis_functions")) {
        basis_set_specs["basis_functions"].get_if(_dofs, "dofs");
      }

      // auto dof_sites_filter = [&_dofs](Site const & _site) ->bool{
      //   if(_dofs.empty() && (_site.dof_size() != 0 || _site.occupant_dof().size() > 1))
      //     return true;
      //
      //   for(DoFKey const &_dof : _dofs) {
      //     if(_site.has_dof(_dof))
      //       return true;
      //     else if(_dof == "occ" && _site.occupant_dof().size() > 1) {
      //       return true;
      //     }
      //   }
      //   return false;
      // };

      // construct cluster orbits  TODO: update with ClusterSpecs
      // make_prim_periodic_orbits(shared_prim, basis_set_specs, dof_sites_filter,
      //                           settings.crystallography_tol(), std::back_inserter(orbits), log());

      // expand the nlist to contain sites in all orbits
      std::set<UnitCellCoord> nbors;
      prim_periodic_neighborhood(orbits.begin(), orbits.end(), std::inserter(nbors, nbors.begin()));
      prim_neighbor_list.expand(nbors.begin(), nbors.end());

      // construct basis functions
      ClexBasis clex_basis {shared_prim, basis_set_specs};
      clex_basis.generate(orbits.begin(), orbits.end(), basis_set_specs);

      // write basis set data files and Clexulator source code
      _write_clexulator(shared_prim, settings, basis_set_name, basis_set_specs, prim_neighbor_list,
                        orbits, clex_basis);

    }
    else {
      // TODO: update
      throw std::runtime_error("Error in `casm bset`: local orbits are not implemented.");
    }
  }

  Clexulator read_clexulator(
    ProjectSettings const &settings,
    std::string const &basis_set_name,
    PrimNeighborList &prim_neighbor_list) {
    return Clexulator {
      settings.project_name() + "_Clexulator",
      settings.dir().clexulator_dir(basis_set_name),
      prim_neighbor_list,
      log(),
      settings.compile_options(),
      settings.so_options()};
  }

}

#define INST_PrimClex_orbits_vec(OrbitType, SymCompareType) \
  template std::back_insert_iterator<std::vector<OrbitType>> PrimClex::orbits( \
    std::string const &, \
    std::back_insert_iterator<std::vector<OrbitType>>, \
    SymCompareType const &) const;

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
