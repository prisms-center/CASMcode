#include <memory>

#include "casm/app/ClexDescription.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/SafeOfstream.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/ChemicalReference.hh"
#include "casm/clex/ClexBasisSpecs.hh"
#include "casm/clex/ClexBasisWriter_impl.hh"
#include "casm/clex/ClexBasis_impl.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/CompositionAxes_impl.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/clex/ECIContainer.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/clex/PrimClex_impl.hh"
#include "casm/clex/io/ProtoFuncsPrinter_impl.hh"
#include "casm/clex/io/file/ChemicalReference_file_io.hh"
#include "casm/clex/io/file/CompositionAxes_file_io.hh"
#include "casm/clex/io/json/ClexBasisSpecs_json_io.hh"
#include "casm/clex/io/json/ECIContainer_json_io.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clusterography/IntegralCluster_impl.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/database/DatabaseHandler_impl.hh"
#include "casm/database/DatabaseTypes_impl.hh"

namespace CASM {

BasicStructure read_prim(ProjectSettings const &project_settings) {
  return read_prim(project_settings.dir().prim(),
                   project_settings.crystallography_tol());
}

std::shared_ptr<Structure const> read_shared_prim(
    ProjectSettings const &project_settings) {
  return std::make_shared<Structure const>(read_prim(project_settings));
}

struct PrimClex::PrimClexData {
  typedef PrimClex::PrimType PrimType;
  typedef std::shared_ptr<PrimType const> PrimType_ptr;

  PrimClexData(ProjectSettings const &_project_settings,
               std::shared_ptr<PrimType const> _shared_prim)
      : settings(_project_settings), prim_ptr(_shared_prim) {
    // Guarantee presence of symmetry info;
    prim_ptr->factor_group();
  }

  PrimClexData(const fs::path &_root)
      : settings(open_project_settings(_root)),
        prim_ptr(read_shared_prim(settings)) {
    // Guarantee presence of symmetry info;
    prim_ptr->factor_group();
  }

  ~PrimClexData() {}

  ProjectSettings settings;

  PrimType_ptr prim_ptr;
  bool vacancy_allowed;
  Index vacancy_index;

  std::unique_ptr<DB::DatabaseHandler> db_handler;

  /// CompositionConverter specifies parameteric composition axes and converts
  /// between
  ///   parametric composition and mol composition
  bool has_composition_axes = false;
  CompositionConverter comp_converter;

  /// ChemicalReference specifies a reference for formation energies, chemical
  /// potentials, etc.
  notstd::cloneable_ptr<ChemicalReference> chem_ref;

  /// Stores the neighboring UnitCell and which sublattices to include in
  /// neighbor lists
  /// - mutable for lazy construction
  mutable std::shared_ptr<PrimNeighborList> nlist;

  typedef std::string BasisSetName;
  mutable std::map<BasisSetName, ClexBasisSpecs> basis_set_specs;
  mutable std::map<BasisSetName, ClexBasis> clex_basis;
  mutable std::map<BasisSetName, Clexulator> clexulator;
  mutable std::map<ClexDescription, ECIContainer> eci;
};

//  **** Constructors ****

/// Initial construction of a PrimClex, from a primitive Structure
PrimClex::PrimClex(ProjectSettings const &_project_settings,
                   std::shared_ptr<PrimType const> _shared_prim)
    : m_data(new PrimClexData(_project_settings, _shared_prim)) {
  m_data->settings.set_crystallography_tol(TOL);

  _init();

  return;
}

/// Construct PrimClex from existing CASM project directory
///  - read PrimClex and directory structure to generate all its Supercells and
///  Configurations, etc.
PrimClex::PrimClex(const fs::path &_root) : m_data(new PrimClexData(_root)) {
  _init();
}

/// Necessary for "pointer to implementation"
PrimClex::~PrimClex() {}

/// Initialization routines
///  - If !root.empty(), read all saved data to generate all Supercells and
///  Configurations, etc.
void PrimClex::_init() {
  auto struc_mol_name = xtal::struc_molecule_name(prim());
  m_data->vacancy_allowed = false;
  for (int i = 0; i < struc_mol_name.size(); ++i) {
    if (xtal::is_vacancy(struc_mol_name[i])) {
      m_data->vacancy_allowed = true;
      m_data->vacancy_index = i;
    }
  }

  if (!has_dir()) {
    return;
  }

  bool read_settings = false;
  bool read_composition = true;
  bool read_chem_ref = true;
  bool read_configs = true;
  bool clear_clex = false;

  refresh(read_settings, read_composition, read_chem_ref, read_configs,
          clear_clex);
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
void PrimClex::refresh(bool read_settings, bool read_composition,
                       bool read_chem_ref, bool read_configs, bool clear_clex) {
  if (read_settings) {
    try {
      m_data->settings = open_project_settings(dir().root_dir());
    } catch (std::exception &e) {
      err_log().error("reading project_settings.json");
      err_log() << "file: " << dir().project_settings() << "\n" << std::endl;
      throw e;
    }
  }

  if (read_composition) {
    m_data->has_composition_axes = false;
    auto comp_axes_path = dir().composition_axes();

    try {
      if (fs::is_regular_file(comp_axes_path)) {
        CompositionAxes opt = read_composition_axes(comp_axes_path);

        if (opt.has_current_axes()) {
          m_data->has_composition_axes = true;
          m_data->comp_converter = opt.curr;
        }
      }
    } catch (std::exception &e) {
      err_log().error("reading composition_axes.json");
      err_log() << "file: " << comp_axes_path << "\n" << std::endl;
      throw e;
    }
  }

  if (read_chem_ref) {
    // read chemical reference
    m_data->chem_ref.reset();
    auto chem_ref_path =
        dir().chemical_reference(m_data->settings.default_clex().calctype,
                                 m_data->settings.default_clex().ref);

    try {
      if (fs::is_regular_file(chem_ref_path)) {
        m_data->chem_ref =
            notstd::make_cloneable<ChemicalReference>(read_chemical_reference(
                chem_ref_path, prim(), settings().lin_alg_tol()));
      }
    } catch (std::exception &e) {
      err_log().error("reading chemical_reference.json");
      err_log() << "file: " << chem_ref_path << "\n" << std::endl;
      throw e;
    }
  }

  if (read_configs) {
    if (m_data->db_handler) {
      // lazy initialization means we just need to close, and the db will be
      // re-opened when needed
      m_data->db_handler->close();
    }
  }

  if (clear_clex) {
    m_data->nlist.reset();
    m_data->clex_basis.clear();
    m_data->clexulator.clear();
    m_data->eci.clear();
  }
}

ProjectSettings &PrimClex::settings() { return m_data->settings; }

const ProjectSettings &PrimClex::settings() const { return m_data->settings; }

bool PrimClex::has_dir() const { return settings().has_dir(); }

const DirectoryStructure &PrimClex::dir() const { return settings().dir(); }

/// \brief Get the crystallography_tol
double PrimClex::crystallography_tol() const { return prim().lattice().tol(); }

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
const PrimClex::PrimType &PrimClex::prim() const { return *(m_data->prim_ptr); }

std::shared_ptr<PrimClex::PrimType const> const &PrimClex::shared_prim() const {
  return this->m_data->prim_ptr;
}

Index PrimClex::n_basis() const { return prim().basis().size(); }

std::shared_ptr<PrimNeighborList> const &PrimClex::shared_nlist() const {
  // lazy neighbor list generation
  if (!m_data->nlist) {
    // construct nlist
    m_data->nlist = std::make_shared<PrimNeighborList>(
        settings().nlist_weight_matrix(),
        settings().nlist_sublat_indices().begin(),
        settings().nlist_sublat_indices().end());
  }

  return m_data->nlist;
}

PrimNeighborList &PrimClex::nlist() const { return *shared_nlist(); }

/// returns true if vacancy are an allowed species
bool PrimClex::vacancy_allowed() const { return m_data->vacancy_allowed; }

/// returns the index of vacancies in composition vectors
Index PrimClex::vacancy_index() const { return m_data->vacancy_index; }

template <typename T>
DB::ValDatabase<T> &PrimClex::generic_db() const {
  return db_handler().template generic_db<T>();
}

template <typename T>
const DB::ValDatabase<T> &PrimClex::const_generic_db() const {
  return db_handler().template const_generic_db<T>();
}

template <typename T>
DB::Database<T> &PrimClex::db() const {
  return db_handler().template db<T>();
}

template <typename T>
const DB::Database<T> &PrimClex::const_db() const {
  return db_handler().template const_db<T>();
}

template <typename T>
DB::PropertiesDatabase &PrimClex::db_props(std::string calc_type) const {
  return db_handler().template db_props<T>(calc_type);
}

template <typename T>
const DB::PropertiesDatabase &PrimClex::const_db_props(
    std::string calc_type) const {
  return db_handler().template const_db_props<T>(calc_type);
}

DB::DatabaseHandler &PrimClex::db_handler() const {
  if (!m_data->db_handler) {
    m_data->db_handler = notstd::make_unique<DB::DatabaseHandler>(*this);
  }
  return *m_data->db_handler;
}

const DB::DatabaseHandler &PrimClex::const_db_handler() const {
  if (!m_data->db_handler) {
    m_data->db_handler = notstd::make_unique<DB::DatabaseHandler>(*this);
  }
  return *m_data->db_handler;
}

bool PrimClex::has_basis_set_specs(std::string const &basis_set_name) const {
  auto it = m_data->basis_set_specs.find(basis_set_name);
  if (it == m_data->basis_set_specs.end()) {
    return fs::exists(dir().bspecs(basis_set_name));
  }
  return true;
}

ClexBasisSpecs const &PrimClex::basis_set_specs(
    std::string const &basis_set_name) const {
  auto it = m_data->basis_set_specs.find(basis_set_name);
  if (it == m_data->basis_set_specs.end()) {
    throw_if_no_basis_set_specs(basis_set_name, dir());

    fs::path basis_set_specs_path = dir().bspecs(basis_set_name);
    jsonParser bspecs_json{basis_set_specs_path};
    ParsingDictionary<DoFType::Traits> const *dof_dict =
        &DoFType::traits_dict();

    InputParser<ClexBasisSpecs> parser{bspecs_json, shared_prim(), dof_dict};
    std::stringstream ss;
    ss << "Error: Invalid file " << basis_set_specs_path;
    report_and_throw_if_invalid(parser, err_log(),
                                std::runtime_error{ss.str()});

    it = m_data->basis_set_specs.emplace(basis_set_name, *parser.value).first;
  }
  return it->second;
}

Clexulator PrimClex::clexulator(std::string const &basis_set_name) const {
  auto it = m_data->clexulator.find(basis_set_name);
  if (it == m_data->clexulator.end()) {
    if (!fs::exists(
            dir().clexulator_src(settings().project_name(), basis_set_name))) {
      std::stringstream ss;
      ss << "Error loading clexulator " << basis_set_name
         << ". No basis functions exist.";
      throw std::runtime_error(ss.str());
    }

    try {
      Clexulator clexulator =
          make_clexulator(settings(), basis_set_name, nlist());
      it = m_data->clexulator.emplace(basis_set_name, clexulator).first;
    } catch (std::exception &e) {
      // TODO: if this fails...
      err_log() << "Error constructing Clexulator. Current settings: \n"
                << std::endl;
      print_compiler_settings_summary(settings(), err_log());

      // TODO: then, try this
      // std::cout << "Error constructing Clexulator. Current settings: \n" <<
      // std::endl; Log tlog(std::cout);
      // print_compiler_settings_summary(settings(), tlog);
      // throw e;
    }
  }
  return it->second;
}

bool PrimClex::has_eci(const ClexDescription &key) const {
  auto it = m_data->eci.find(key);
  if (it == m_data->eci.end()) {
    return fs::exists(
        dir().eci(key.property, key.calctype, key.ref, key.bset, key.eci));
  }
  return true;
}

const ECIContainer &PrimClex::eci(const ClexDescription &key) const {
  auto it = m_data->eci.find(key);
  if (it == m_data->eci.end()) {
    fs::path eci_path =
        dir().eci(key.property, key.calctype, key.ref, key.bset, key.eci);
    if (!fs::exists(eci_path)) {
      std::stringstream ss;
      ss << "Error loading ECI. eci.json does not exist.\n  Expected at: "
         << eci_path.string();
      throw std::runtime_error(ss.str());
    }

    jsonParser eci_json{eci_path};
    InputParser<ECIContainer> parser{eci_json};
    std::stringstream ss;
    ss << "Error loading ECI: Invalid file " << eci_path;
    report_and_throw_if_invalid(parser, err_log(),
                                std::runtime_error{ss.str()});

    it = m_data->eci.emplace(key, *parser.value).first;
  }
  return it->second;
}

/// Write clust.json, basis.json, and clexulator source code, given orbits
struct WriteBasisSetDataImpl {
  WriteBasisSetDataImpl(std::shared_ptr<Structure const> _shared_prim,
                        ProjectSettings const &_settings,
                        std::string const &_basis_set_name,
                        ClexBasisSpecs const &_basis_set_specs,
                        PrimNeighborList &_prim_neighbor_list)
      : shared_prim(_shared_prim),
        settings(_settings),
        basis_set_name(_basis_set_name),
        basis_set_specs(_basis_set_specs),
        prim_neighbor_list(_prim_neighbor_list) {}

  std::shared_ptr<Structure const> shared_prim;
  ProjectSettings const &settings;
  std::string const &basis_set_name;
  ClexBasisSpecs const &basis_set_specs;
  PrimNeighborList &prim_neighbor_list;

  template <typename OrbitVecType>
  void operator()(OrbitVecType const &orbits) const {
    const auto &dir = settings.dir();
    ParsingDictionary<DoFType::Traits> const *dof_dict =
        &DoFType::traits_dict();
    double xtal_tol = settings.crystallography_tol();

    // generate ClexBasis
    ClexBasis clex_basis{shared_prim, basis_set_specs, dof_dict};
    clex_basis.generate(orbits.begin(), orbits.end());

    // delete any existing data
    throw_if_no_basis_set_specs(basis_set_name, dir);
    dir.delete_bset_data(settings.project_name(), basis_set_name);

    jsonParser basis_set_specs_json;
    to_json(basis_set_specs, basis_set_specs_json, *shared_prim, dof_dict);

    // write clust
    fs::path clust_json_path = dir.clust(basis_set_name);
    jsonParser clust_json;
    OrbitPrinterOptions orbit_printer_options;
    orbit_printer_options.print_invariant_group = true;
    ProtoSitesPrinter sites_printer{orbit_printer_options};
    write_clust(orbits.begin(), orbits.end(), clust_json, sites_printer,
                basis_set_specs_json);
    clust_json.write(clust_json_path);

    // write basis
    fs::path basis_json_path = dir.basis(basis_set_name);
    jsonParser basis_json;
    write_site_basis_funcs(shared_prim, clex_basis, basis_json);
    bool align = false;
    ProtoFuncsPrinter funcs_printer{clex_basis, shared_prim->shared_structure(),
                                    align, orbit_printer_options};
    write_clust(orbits.begin(), orbits.end(), basis_json, funcs_printer,
                basis_set_specs_json);
    basis_json.write(basis_json_path);

    // write source code
    fs::path clexulator_src_path =
        dir.clexulator_src(settings.project_name(), basis_set_name);
    std::string clexulator_name = settings.global_clexulator_name();
    auto param_pack_type = basis_set_specs.basis_function_specs.param_pack_type;
    fs::ofstream outfile;
    outfile.open(clexulator_src_path);
    ClexBasisWriter clexwriter{*shared_prim, param_pack_type};
    clexwriter.print_clexulator(clexulator_name, clex_basis, orbits,
                                prim_neighbor_list, outfile, xtal_tol);
    outfile.close();
  }
};

/// Write clexulator, clust.json, basis.json by constructing orbits and
/// ClexBasis
void write_basis_set_data(std::shared_ptr<Structure const> shared_prim,
                          ProjectSettings const &settings,
                          std::string const &basis_set_name,
                          ClexBasisSpecs const &basis_set_specs,
                          PrimNeighborList &prim_neighbor_list) {
  auto const &cluster_specs = *basis_set_specs.cluster_specs;
  Log &log = CASM::log();

  WriteBasisSetDataImpl writer{shared_prim, settings, basis_set_name,
                               basis_set_specs, prim_neighbor_list};
  for_all_orbits(cluster_specs, log, writer);
}

/// Make Clexulator from existing source code
///
/// Notes:
/// - Use `write_basis_set_data` to write Clexulator source code prior to
/// calling this function.
/// - This function will compile the Clexulator source code if that has not yet
/// been done.
Clexulator make_clexulator(ProjectSettings const &settings,
                           std::string const &basis_set_name,
                           PrimNeighborList &prim_neighbor_list) {
  throw_if_no_clexulator_src(settings.project_name(), basis_set_name,
                             settings.dir());
  return Clexulator{settings.project_name() + "_Clexulator",
                    settings.dir().clexulator_dir(basis_set_name),
                    prim_neighbor_list, settings.compile_options(),
                    settings.so_options() + " -lcasm "};
}

}  // namespace CASM

#include "casm/database/DatabaseTypes.hh"

// explicit template instantiations
#define INST_PrimClex(r, data, type)                                       \
  template DB::Database<type> &PrimClex::db<type>() const;                 \
  template const DB::Database<type> &PrimClex::const_db<type>() const;     \
  template DB::ValDatabase<type> &PrimClex::generic_db<type>() const;      \
  template const DB::ValDatabase<type> &PrimClex::const_generic_db<type>() \
      const;

// explicit template instantiations
#define INST_PrimClexProps(r, data, type)                                \
  template DB::PropertiesDatabase &PrimClex::db_props<type>(             \
      std::string calc_type) const;                                      \
  template const DB::PropertiesDatabase &PrimClex::const_db_props<type>( \
      std::string calc_type) const;

namespace CASM {
BOOST_PP_SEQ_FOR_EACH(INST_PrimClex, _, CASM_DB_TYPES)
BOOST_PP_SEQ_FOR_EACH(INST_PrimClexProps, _, CASM_DB_CONFIG_TYPES)
}  // namespace CASM
