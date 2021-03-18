#ifndef CASM_Configuration
#define CASM_Configuration

#include <map>
#include <memory>

#include "casm/clex/Calculable.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/ConfigurationTraits.hh"
#include "casm/clex/HasCanonicalForm.hh"
#include "casm/clex/HasSupercell.hh"
#include "casm/clusterography/ClusterDecl.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/database/Cache.hh"
#include "casm/database/Database.hh"
#include "casm/database/Named.hh"
#include "casm/misc/Comparisons.hh"
#include "casm/symmetry/PermuteIterator.hh"

namespace CASM {
namespace xtal {
class Molecule;
class UnitCellCoord;
}  // namespace xtal
using xtal::Molecule;
using xtal::UnitCellCoord;

class IntegralCluster;
class PrimClex;
class Supercell;
class Clexulator;
class ConfigIsEquivalent;
template <typename ConfigType, typename IsEqualImpl>
class GenericConfigCompare;
using ConfigCompare = GenericConfigCompare<Configuration, ConfigIsEquivalent>;

struct ConfigInsertResult;
struct RefToCanonicalPrim;

/// \defgroup Configuration
///
/// \brief A Configuration represents the values of all degrees of freedom in a
/// Supercell
///
/// \ingroup Clex
///
/// @{

typedef ConfigCanonicalForm<
    HasSupercell<Comparisons<Calculable<CRTPBase<Configuration>>>>>
    ConfigurationBase;

/// Configuration, a periodic perturbation of the infinite crystal within the
/// DoF space of the prim
///
/// A Configuration represents a particular periodic perturbation of the
/// infinite crystal within
///   the space of allowed perturbations defined by the BasicStructure.
///   Configuration has:
/// - a Supercell, representing the translational periodicity of the
/// perturbation
/// - a ConfigDoF, representing the values of DoF (local discrete and local
/// continuous) "within"
//    the supercell (i.e. the translationally unique perturbation), plus the
//    values of the
///   global DoF.
/// - calculated properties (std::map<std::string, MappedProperties>), the
/// values of properties
///   (i.e. energy, relaxation displacements, relaxation strain, etc.) that are
///   dependent on the values of the DoF. Properties are stored in a map with
///   key "calctype" to allow for different values of the properties depending
///   on the calculation method.
///
class Configuration : public ConfigurationBase {
 public:
  //********* CONSTRUCTORS *********

  /// Construct a default Configuration, with a shared Supercell
  ///
  /// Note:
  /// - This is one of the preferred constructors for new code. It can be used
  /// when the Supercell
  ///   is not stored in a Database<Supercell>. In the future, it will be used
  ///   in that context also.
  /// - This Configuration does own its own Supercell, so the lifetime of the
  /// Supercell is
  ///   guaranteed to exceed the lifetime of this Configuration
  explicit Configuration(
      std::shared_ptr<Supercell const> const &_supercell_ptr);

  /// Construct a default Configuration, with a shared Supercell
  ///
  /// Note:
  /// - This is one of the preferred constructors for new code. It can be used
  /// when the Supercell
  ///   is not stored in a Database<Supercell>. In the future, it will be used
  ///   in that context also.
  /// - This Configuration does own its own Supercell, so the lifetime of the
  /// Supercell is
  ///   guaranteed to exceed the lifetime of this Configuration
  explicit Configuration(std::shared_ptr<Supercell const> const &_supercell_ptr,
                         ConfigDoF const &_dof);

  /// Build a Configuration sized to _supercell with all fields initialized and
  /// set to zero
  static Configuration zeros(
      const std::shared_ptr<Supercell const> &_supercell_ptr);

  /// Build a Configuration sized to _supercell with all fields initialized and
  /// set to zero
  static Configuration zeros(
      const std::shared_ptr<Supercell const> &_supercell_ptr, double _tol);

  // *** The following constructors should be avoided in new code, if possible
  //
  //     They are currently still required in some code. In the future:
  //     - Configuration will make exclusive use of std::shared_ptr<Supercell
  //     const>
  //     - The jsonConstructor<Configuration>::from_json method will be used to
  //       construct Configuration from JSON

  /// Construct a default Configuration, with a pointer to a Supercell and
  ///
  /// Note:
  /// - Whenever possible, this constructor should not be used in new code .
  /// - This constructor keeps a pointer to _supercell, whose lifetime must
  /// exceed the lifetime
  ///   of this Configuration
  explicit Configuration(Supercell const &_supercell);

  /// Construct a default Configuration
  ///
  /// Note:
  /// - Whenever possible, this constructor should not be used in new code .
  /// - This constructor keeps a pointer to _supercell, whose lifetime must
  /// exceed the lifetime
  ///   of this Configuration
  explicit Configuration(const Supercell &_supercell, const ConfigDoF &_dof);

  /// Build a Configuration sized to _scel with all fields initialized and set
  /// to zero
  ///
  /// Note:
  /// - Whenever possible, this constructor should not be used in new code .
  static Configuration zeros(Supercell const &_scel);

  /// Build a Configuration sized to _scel with all fields initialized and set
  /// to zero
  ///
  /// Note:
  /// - Whenever possible, this constructor should not be used in new code .
  static Configuration zeros(Supercell const &_scel, double _tol);

  // ******** Supercell **********************

  /// \brief Get the Supercell for this Configuration
  const Supercell &supercell() const;

  const Lattice &ideal_lattice() const;

  /// Returns number of sites, NOT the number of primitives that fit in here
  Index size() const;

  /// \brief Get the UnitCellCoord for a given linear site index
  UnitCellCoord uccoord(Index site_l) const;

  /// \brief Return the linear index corresponding to integral coordinates
  Index linear_index(const UnitCellCoord &bijk) const;

  /// \brief Get the basis site index for a given linear linear site index
  int sublat(Index site_l) const;

  // ******** Degrees of Freedom **************
  //
  // ** Note: Calculated properties are not automatically updated when dof are
  // changed, **
  // **       nor are the written records automatically updated **

  /// \brief const Access the DoF
  const ConfigDoF &configdof() const { return m_configdof; }

  /// \brief Access the DoF
  ///
  /// - This will invalidate the Configuration's id
  ConfigDoF &configdof() {
    _modify_dof();
    return m_configdof;
  }

  /// \brief Clear all DoF
  ///
  /// - This will invalidate the Configuration's id
  void clear();

  // ----- Occupation ------------

  /// \brief Set occupant variables to background structure
  ///
  /// - This will invalidate the Configuration's id
  void init_occupation();

  /// \brief Set occupant variables
  ///
  /// With one value for each site in the Configuration, this std::vector
  /// describes which occupant is at each of the 'N' sites of the configuration.
  /// The occupant on site l can be obtained from the occupation variable using:
  /// \code
  /// Molecule on_site_l = config.prim().basis()[ config.sublat(l)
  /// ].site_occupant[ config.occupation()[l]]; \endcode
  /// - For a CASM project, the occupation variables will be ordered according
  /// to the occupant DoF in a "prim.json" file. This means that for the
  /// background structure, 'occupation' is all 0
  /// - This will invalidate the Configuration's id
  ///
  /// \throws If \code newoccupation.size() != this->size() \endcode
  ///
  /// set_occupation ensures that ConfigDoF::size() is compatible with
  /// _occupation.size() or if ConfigDoF::size()==0, sets ConfigDoF::size() to
  /// _occupation.size()
  void set_occupation(Eigen::Ref<const Eigen::VectorXi> const &_occupation) {
    configdof().set_occupation(_occupation);
  }

  /// \brief Occupant variables
  ///
  /// With one value for each site in the Configuration, this std::vector
  /// describes which occupant is at each of the 'N' sites of the configuration.
  /// The occupant on site l can be obtained from the occupation variable using:
  /// \code
  /// Molecule on_site_l = config.prim().basis()[ config.sublat(l)
  /// ].site_occupant[ config.occupation()[l]]; \endcode
  /// - For a CASM project, the occupation variables will be ordered according
  /// to the occupant DoF in a "prim.json" file. This means that for the
  /// background structure, 'occupation' is all 0
  ///
  Eigen::VectorXi const &occupation() const { return configdof().occupation(); }

  /// \brief Set occupant variable on site l
  ///
  /// The occupant on site l can be obtained from the occupation variable using:
  /// \code
  /// Molecule on_site_l = config.prim().basis()[ config.sublat(l)
  /// ].site_occupant[config.occ(l)]; \endcode
  /// - For a CASM project, the occupation variables will be ordered according
  /// to the occupant DoF in a "prim.json" file. This means that for the
  /// background structure, 'occupation' is all 0
  /// - No check is performed as to whether val is in the correct range
  /// - This will invalidate the Configuration's id
  ///
  void set_occ(Index site_l, int val);

  /// \brief Occupant variable on site l
  ///
  /// The occupant on site l can be obtained from the occupation variable using:
  /// \code
  /// Molecule on_site_l = config.prim().basis()[ config.sublat(l)
  /// ].site_occupant[config.occ(l)]; \endcode
  /// - For a CASM project, the occupation variables will be ordered according
  /// to the occupant DoF in a "prim.json" file. This means that for the
  /// background structure, 'occupation' is all 0
  ///
  const int &occ(Index site_l) const { return configdof().occ(site_l); }

  /// \brief Molecule on site l
  ///
  /// Equivalent to:
  /// \code
  /// config.prim().basis()[ config.sublat(l) ].site_occupant[
  /// config.occupation()[l]]; \endcode
  ///
  const Molecule &mol(Index site_l) const;

  /// \brief True if Configuration has occupation DoF
  bool has_occupation() const { return configdof().has_occupation(); }

  /// \brief Clear occupation
  ///
  /// - This will invalidate the Configuration's id
  // void clear_occupation();

  // ******** Comparisons, Symmetry, Crystallography  ******

  /// \brief Compare Configuration, via ConfigCompare
  ///
  /// - Must have the same Supercell
  bool operator<(const Configuration &B) const;

  ConfigCompare less() const;

  ConfigIsEquivalent equal_to() const;

  /// \brief Check if this is a primitive Configuration
  bool is_primitive() const;

  /// \brief Returns a PermuteIterator corresponding to the first non-zero pure
  /// translation that maps the Configuration onto itself.
  PermuteIterator find_translation() const;

  /// \brief Return the primitive Configuration
  Configuration primitive() const;

  /// \brief Returns the canonical form Configuration in the canonical Supercell
  Configuration in_canonical_supercell() const;

  /// \brief Returns the subgroup of the Supercell factor group that leaves the
  ///        Configuration unchanged
  std::vector<PermuteIterator> factor_group() const;

  /// \brief Gives the subgroup of the supercell that leaves this configuration
  /// unchanged
  using ConfigurationBase::invariant_subgroup;
  std::vector<PermuteIterator> invariant_subgroup() const;

  /// \brief Determines if this Configuration is in canonical form
  using ConfigurationBase::is_canonical;
  bool is_canonical() const;

  /// \brief Get symmetric multiplicity, excluding translations
  int multiplicity() const;

  /// \brief Returns the point group that leaves the Configuration unchanged
  std::vector<PermuteIterator> point_group() const;

  /// \brief Returns the point group that leaves the Configuration unchanged
  std::string point_group_name() const;

  /// \brief Transform Configuration from PermuteIterator via *this =
  /// permute_iterator * *this
  Configuration &apply_sym(const PermuteIterator &it);

  //********** Composition ***********

  // Returns composition on each sublattice: sublat_comp[ prim basis site /
  // sublattice][ molecule_type]
  //   molucule_type is ordered as in the Prim structure's site_occupant list
  //   for that basis site (includes vacancies)
  std::vector<Eigen::VectorXd> sublattice_composition() const;

  // Returns number of each molecule by sublattice:
  //   sublat_num_each_molecule[ prim basis site / sublattice ][ molecule_type]
  //   molucule_type is ordered as in the Prim structure's site_occupant list
  //   for that basis site
  std::vector<Eigen::VectorXi> sublat_num_each_molecule() const;

  // Returns composition, not counting vacancies
  //    composition[ molecule_type ]: molecule_type ordered as prim structure's
  //    get_struc_molecule(), with [Va]=0.0
  Eigen::VectorXd composition() const;

  // Returns composition, including vacancies
  //    composition[ molecule_type ]: molecule_type ordered as prim structure's
  //    get_struc_molecule()
  Eigen::VectorXd true_composition() const;

  /// Returns num_each_molecule[ molecule_type], where 'molecule_type' is
  /// ordered as Structure::get_struc_molecule()
  Eigen::VectorXi num_each_molecule() const;

  /// Returns parametric composition, as calculated using PrimClex::param_comp
  Eigen::VectorXd param_composition() const;

  /// Returns num_each_component[ component_type] per prim cell,
  ///   where 'component_type' is ordered as ParamComposition::components
  Eigen::VectorXd num_each_component() const;

  //********* IO ************

  /// \brief Insert this configuration (in primitive & canonical form) in the
  /// database
  ConfigInsertResult insert(bool primitive_only = false) const;

  // /// Writes the Configuration to JSON
  // jsonParser &to_json(jsonParser &json) const;
  //
  // /// Reads the Configuration from JSON
  // void from_json(const jsonParser &json, const Supercell &scel,
  //                std::string _id);
  //
  // /// Reads the Configuration from JSON
  // void from_json(const jsonParser &json, const PrimClex &primclex,
  //                std::string _configname);

  /// \brief Split configuration name string into scelname and config id
  static std::pair<std::string, std::string> split_name(std::string configname);

 private:
  friend struct Comparisons<Calculable<CRTPBase<Configuration>>>;
  friend class DB::Named<CRTPBase<Configuration>>;

  std::string generate_name_impl() const;

  /// \brief operator== comparison of Configuration, via ConfigEqual
  ///
  /// - Must have the same Supercell
  /// - Checks that all DoF are the same, within tolerance
  bool eq_impl(const Configuration &B) const;

  /// Pointer to the Supercell for this Configuration
  Supercell const *m_supercell;

  /// Used when constructing temporary Configuration in non-canonical Supercell
  std::shared_ptr<Supercell const> m_supercell_ptr;

  /// Degrees of Freedom
  ConfigDoF m_configdof;
};

/// \brief Holds results of Configuration::insert
///
/// - 'canonical' refers to the canonical form of the Configuration in it's
///   canonical equivalent Supercell.  The canonical form may be primitive or
///   non-primitive
/// - 'primitive' refers to the primitive canonical Configuration.
///
struct ConfigInsertResult {
  typedef DB::DatabaseIterator<Configuration> iterator;

  /// True if primitive did not exist before insertion
  bool insert_primitive;

  /// Iterator pointing at primitive
  iterator primitive_it;

  /// True if canonical configuration did not exist before insertion
  bool insert_canonical;

  /// Iterator pointing at canonical, if existing
  iterator canonical_it;
};

/// \brief Operations that transform a canonical primitive configuration to any
/// equivalent
///
/// \code
/// Configuration equiv_prim_config = copy_apply(from_canonical_config,
/// prim_canon_config); FillSupercell f {config.supercell()}; Configuration
/// config = f(from_canonical_lat, equiv_prim_config); \endcode
struct RefToCanonicalPrim {
  /// \brief Get operations that transform canonical primitive to this
  RefToCanonicalPrim(const Configuration &_config);

  std::string name() const;

  Configuration config;
  Configuration prim_canon_config;
  SymOp from_canonical_lat;
  PermuteIterator from_canonical_config;
  Eigen::Matrix3i transf_mat;
};

/// \brief Returns the sub-configuration that fills a particular Supercell
///
/// \param sub_scel_ptr The Supercell of the sub-configuration
/// \param super_config The super-configuration
/// \param origin The UnitCell indicating the which unit cell in the
///        super-configuration is the origin in sub-configuration
///
/// - Copies DoF from the super-configuration directly into the
/// sub-configuration
///
Configuration sub_configuration(std::shared_ptr<Supercell const> sub_scel_ptr,
                                const Configuration &super_config,
                                const UnitCell &origin = UnitCell(0, 0, 0));

/// \brief Make Configuration from name string
Configuration make_configuration(const PrimClex &primclex, std::string name);

/// \brief Returns correlations using 'clexulator'.
Eigen::VectorXd correlations(const Configuration &config,
                             Clexulator const &clexulator);

/// Returns correlation contribution from a single unit cell, not normalized.
Eigen::VectorXd corr_contribution(Index linear_unitcell_index,
                                  const Configuration &config,
                                  Clexulator const &clexulator);

/// \brief Returns point correlations from a single site, normalized by cluster
/// orbit size
Eigen::VectorXd point_corr(Index linear_unitcell_index, Index neighbor_index,
                           const Configuration &config,
                           Clexulator const &clexulator);

/// \brief Returns gradient correlations using 'clexulator', with respect to DoF
/// 'dof_type'
Eigen::MatrixXd gradcorrelations(const Configuration &config,
                                 Clexulator const &clexulator, DoFKey &key);

/// Returns parametric composition, as calculated using PrimClex::param_comp
Eigen::VectorXd comp(const Configuration &config);

/// \brief Returns the composition, as number of each species per unit cell
Eigen::VectorXd comp_n(const Configuration &config);

/// \brief Returns the vacancy composition, as number per unit cell
double n_vacancy(const Configuration &config);

/// \brief Returns the total number species per unit cell
double n_species(const Configuration &config);

/// \brief Returns the composition as species fraction, with [Va] = 0.0, in
///        the order of Structure::get_struc_molecule
Eigen::VectorXd species_frac(const Configuration &config);

/// \brief Returns the composition as site fraction, in the order of
/// Structure::get_struc_molecule
Eigen::VectorXd site_frac(const Configuration &config);

/// \brief Returns the energy, normalized per unit cell
double energy(const Configuration &config);

/// \brief Returns the energy, normalized per species
double energy_per_species(const Configuration &config);

/// \brief Returns the reference energy, normalized per unit cell
double reference_energy(const Configuration &config);

/// \brief Returns the reference energy, normalized per species
double reference_energy_per_species(const Configuration &config);

/// \brief Returns the formation energy, normalized per unit cell
double formation_energy(const Configuration &config);

/// \brief Returns the formation energy, normalized per species
double formation_energy_per_species(const Configuration &config);

/// \brief Returns the formation energy, normalized per unit cell
double clex_formation_energy(const Configuration &config);

/// \brief Returns the formation energy, normalized per species
double clex_formation_energy_per_species(const Configuration &config);

/// \brief Root-mean-square forces of relaxed configurations, determined from
/// DFT (eV/Angstr.)
double rms_force(const Configuration &_config);

/// \brief Cost function that describes the degree to which basis sites have
/// relaxed
double atomic_deformation(const Configuration &_config);

/// \brief Cost function that describes the degree to which lattice has relaxed
double lattice_deformation(const Configuration &_config);

/// \brief Change in volume due to relaxation, expressed as the ratio V/V_0
double volume_relaxation(const Configuration &_config);

/// \brief Returns the relaxed magnetic moment, normalized per unit cell
double relaxed_magmom(const Configuration &_config);

/// \brief Returns the relaxed magnetic moment, normalized per species
double relaxed_magmom_per_species(const Configuration &_config);

/// \brief relaxed forces of configuration, determined from DFT (eV/Angstr.), as
/// a 3xN matrix
Eigen::MatrixXd relaxed_forces(const Configuration &_config);

/// \brief relaxed forces of configuration, determined from DFT (eV/Angstr.), as
/// a 3xN matrix
Eigen::MatrixXd relaxed_forces(const Configuration &_config);

/// \brief Returns an Integral Cluster representing the perturbed sites between
/// the configs
IntegralCluster config_diff(const Configuration &_config1,
                            const Configuration &_config2);

/// Returns a rotated/translated version of config 2 that leaves it closest to
/// the occupation of config1
Configuration closest_setting(const Configuration &_config1,
                              const Configuration &_config2);

/// \brief Returns a Configuration with the sites in _clust clipped from _config
/// and placed in _bg
Configuration config_clip(const Configuration &_config,
                          const Configuration &_bg, IntegralCluster &_clust);

/// \brief returns true if _config describes primitive cell of the configuration
/// it describes
bool is_primitive(const Configuration &_config);

/// \brief returns true if no symmetry transformation applied to _config will
/// increase its lexicographic order
bool is_canonical(const Configuration &_config);

/// \brief returns true if _config is an endpoint of any existing
/// diff_trans_config
bool is_diff_trans_endpoint(const Configuration &_config);

/// \brief returns which diff_trans _config is an endpoint of
std::string diff_trans_endpoint_of(const Configuration &_config);

bool has_energy(const Configuration &_config);

bool has_reference_energy(const Configuration &_config);

bool has_formation_energy(const Configuration &_config);

bool has_rms_force(const Configuration &_config);

bool has_atomic_deformation(const Configuration &_config);

bool has_lattice_deformation(const Configuration &_config);

bool has_volume_relaxation(const Configuration &_config);

bool has_relaxed_magmom(const Configuration &_config);

bool has_relaxed_mag_basis(const Configuration &_config);

/// \brief Returns correlations using 'clexulator'. Supercell needs a correctly
/// populated neighbor list.
Eigen::VectorXd correlations(const ConfigDoF &configdof, const Supercell &scel,
                             Clexulator const &clexulator);

/// Returns correlation contribution from a single unit cell, not normalized.
Eigen::VectorXd corr_contribution(Index linear_unitcell_index,
                                  const ConfigDoF &configdof,
                                  const Supercell &scel,
                                  Clexulator const &clexulator);

/// \brief Returns point correlations from a single site, normalized by cluster
/// orbit size
Eigen::VectorXd point_corr(Index linear_unitcell_index, Index neighbor_index,
                           const ConfigDoF &configdof, const Supercell &scel,
                           Clexulator const &clexulator);

/// \brief Returns gradient correlations using 'clexulator', with respect to DoF
/// 'dof_type'
Eigen::MatrixXd gradcorrelations(const ConfigDoF &configdof,
                                 const Supercell &scel,
                                 Clexulator const &clexulator, DoFKey &key);

/// \brief Returns num_each_molecule(molecule_type), where 'molecule_type' is
/// ordered as Structure::get_struc_molecule()
Eigen::VectorXi num_each_molecule(const ConfigDoF &configdof,
                                  const Supercell &scel);

/// \brief Returns comp_n, the number of each molecule per primitive cell,
/// ordered as Structure::get_struc_molecule()
Eigen::VectorXd comp_n(const ConfigDoF &configdof, const Supercell &scel);

/** @} */

}  // namespace CASM

#endif
