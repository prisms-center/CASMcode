#ifndef CASM_ConfigMapping
#define CASM_ConfigMapping

#include <vector>
#include <unordered_set>
#include "casm/CASM_global_definitions.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
struct HermiteHash {
  std::size_t operator()(Eigen::Matrix3i mat) const {
    std::string unroll_str;
    Eigen::VectorXi vec = CASM::HermiteCounter_impl::_canonical_unroll(mat);
    for(int i = 0; i < vec.size(); ++i) {
      unroll_str = unroll_str + std::to_string(vec(i)) + ",";
    }
    std::size_t const h(std::hash<std::string> {}(unroll_str));
    return h;
  }
};
namespace CASM {
  class Supercell;
  class Lattice;
  class SymGroup;
  class PermuteIterator;
  class PrimClex;
  class Configuration;
  class ConfigDoF;
  class Structure;

  namespace Completer {
    class ImportOption;
  }

  /// Data structure holding results of ConfigMapper algorithm
  struct ConfigMapperResult {

    ConfigMapperResult() :
      success(false) {}

    /// Output structure, after applying lattice similarity and/or rotation to
    /// input structure.
    SimpleStructure structure;

    /// The configuration the input structure was mapped onto
    std::unique_ptr<Configuration> config;

    /// relaxation_properties is populated by relaxation properties:
    ///
    /// - 'lattice_deformation': lattice mapping score
    /// - 'basis_deformation': atomic mapping score
    /// - 'volume_relaxation': V/V_ideal
    /// - 'relaxation_deformation': 3x3 tensor describing cell relaxation
    /// - 'relaxation_displacement': Nx3 matrix describing basis displacements
    /// - 'relaxed_energy': the energy of the relaxed configuration
    jsonParser relaxation_properties;

    /// best_assignment is populated by the permutation of sites in the imported
    /// structure that maps them onto sites of the ideal crystal (excluding vacancies)
    std::vector<Index> best_assignment;

    /// cart_op is populated by the cartesian isometry that rotates the imported
    /// structure onto the coordinate system of the ideal crystal
    Eigen::Matrix3d cart_op;

    /// Rigid translation isometry that puts the imported structures origin onto the origin of ideal crystal
    Eigen::Vector3d trans;

    /// True if could map to prim, false if not
    bool success;

    /// Failure message if could not map to prim
    std::string fail_msg;

  };

  /// Data structure holding a potential mapping, which consists of deformation, occupation array, and displacement field
  struct MappedConfig {
    // typedefs to provide flexibility if we eventually change to a Eigen::Matrix<3,Eigen::Dynamic>
    typedef Eigen::MatrixXd DisplacementMatrix;

    // Can treat as a Eigen::VectorXd
    typedef DisplacementMatrix::ColXpr Displacement;
    typedef DisplacementMatrix::ConstColXpr ConstDisplacement;

    MappedConfig() {
      deformation.setIdentity();
    }

    MappedConfig(Index Nsites) :
      deformation(Eigen::Matrix3d::Identity()),
      occupation(Eigen::VectorXi::Zero(Nsites)),
      displacement(Eigen::MatrixXd::Zero(3, Nsites)) {}

    void clear() {
      throw std::runtime_error("MappedConfig::clear() not implemented");
    }

    Displacement disp(Index i) {
      return displacement.col(i);
    }

    ConstDisplacement disp(Index i) const {
      return displacement.col(i);
    }

    ConfigDoF to_configdof() const;

    MappedConfig &apply_sym(PermuteIterator const &it);


    Eigen::Matrix3d deformation;
    Eigen::VectorXi occupation;
    Eigen::MatrixXd displacement;
  };

  /// A class for mapping an arbitrary crystal structure as a configuration of a crystal template
  /// as described by a PrimClex.  ConfigMapper manages options for the mapping algorithm and mapping cost function
  /// It also caches some information about supercell lattices so that batch imports are more efficient
  ///
  /// \ingroup Configuration
  class ConfigMapper {
  public:
    typedef std::map<Index, std::vector<Lattice> > LatMapType;
    enum NullInitializer {null_initializer};
    enum Options {none = 0,
                  rotate = (1u << 0),
                  strict = (1u << 1),
                  robust = (1u << 2)
                 };

    ///\brief Default construction not allowed -- this constructor provides an override
    ConfigMapper(NullInitializer) :
      m_pclex(nullptr),
      m_lattice_weight(0.5),
      m_max_volume_change(0.5),
      m_min_va_frac(0.),
      m_max_va_frac(1.) {
    }

    ///\brief Construct and initialize a ConfigMapper
    ///\param _pclex the PrimClex that describes the crystal template
    ///
    ///\param _lattice_weight
    ///\parblock
    ///          free parameter 'w' in the cost function: total_cost = w*lattice_deformation+(1-w)*basis_deformation
    ///          can vary between 0 (completely basis-focused) and 1 (completely lattice-focused)
    ///\endparblock
    ///
    ///\param _max_volume_change
    ///\parblock
    ///          constrains the search space by assuming a limit on allowed volume change
    ///          only taken into account when non-interstitial vacancies are allowed
    ///\endparblock
    ///
    ///\param _options
    ///\parblock
    ///          specify a combination of ConfigMapper::Options using bitwise OR: Ex. _options=ConfigMapper::rotate|ConfigMapper::strict
    ///          Options are:
    ///             'rotate': removes rigid rotation of the imported crystal, in a least-squares sense (i.e., yields a symmetric deformation tensor)
    ///             'robust': does not assume the imported structure might be ideal ('robust' is much slower for importing ideal structures, but if 'robust' is not
    ///                       set and a non-ideal structure is passed, this will be almost always be detected and robust methods will be used instead. Thus, 'robust'
    ///                       is slightly faster if imported Structures are *not* ideal)
    ///             'strict': prevents transformation into canonical form. Tries to preserve original orientation of imported structure if possible
    ///\endparblock
    ///
    ///\param _tol tolerance for mapping comparisons
    ConfigMapper(const PrimClex &_pclex,
                 double _lattice_weight,
                 double _max_volume_change = 0.5,
                 int _options = robust, // this should actually be a bitwise-OR of ConfigMapper::Options
                 double _tol = TOL);


    const PrimClex &primclex() const {
      return *m_pclex;
    }

    void set_primclex(const PrimClex &_pclex) {
      m_pclex = &_pclex;
    }

    double lattice_weight() const {
      return m_lattice_weight;
    }

    void set_lattice_weight(double _lw) {
      m_lattice_weight = max(min(_lw, 1.0), 1e-9);
    }

    double min_va_frac() const {
      return m_min_va_frac;
    }

    void set_min_va_frac(double _min_va) {
      m_min_va_frac = max(_min_va, 0.);
    }

    double max_va_frac() const {
      return m_max_va_frac;
    }

    void set_max_va_frac(double _max_va) {
      m_max_va_frac = min(_max_va, 1.);
    }


    ///\brief imports structure specified by 'pos_path' into primclex() by finding optimal mapping
    ///       and then setting displacements and strain to zero (only the mapped occupation is preserved)
    ///
    ConfigMapperResult import_structure_occupation(const fs::path &pos_path) const;

    ///\brief imports structure specified by '_struc' into primclex()
    ///
    ConfigMapperResult import_structure_occupation(const SimpleStructure &_struc) const;

    ///\brief imports structure specified by '_struc' into primclex()
    ///\param hint_ptr[in]
    ///\parblock
    ///                provides a suggestion for which Configuration _struc should map onto
    ///                The hint is used to reduce search times, but may be used in the future
    ///                in combination with Option 'strict' to force mapping onto a particular configuration
    ///                or be used to provide user reports of the form "Suggested mapping: 0.372; Optimal mapping: 0.002"
    ///\endparblock
    ///
    ConfigMapperResult import_structure_occupation(const SimpleStructure &_struc,
                                                   const Configuration *hint_ptr) const;


    ///\brief imports structure specified by 'pos_path' into primclex() by finding optimal mapping
    ///       unlike import_structure_occupation, displacements and strain are preserved
    ///
    ConfigMapperResult import_structure(const fs::path &pos_path) const;

    ///\brief imports structure specified by '_struc' into primclex() by finding optimal mapping
    ///       unlike import_structure_occupation, displacements and strain are preserved
    ///
    ConfigMapperResult import_structure(const SimpleStructure &_struc) const;


    ///\brief Low-level routine to map a structure onto a MappedConfig
    ///\param mapped_config[out] MappedConfig that is result of mapping procedure
    ///\param mapped_lat[out] Ideal supercell lattice (in Niggli form) of mapped configuration
    bool struc_to_configdof(const SimpleStructure &_struc,
                            MappedConfig &mapped_config,
                            Lattice &mapped_lat) const;

    ///\brief Low-level routine to map a structure onto a ConfigDof
    ///\param mapped_config[out] MappedConfig that is result of mapping procedure
    ///\param mapped_lat[out] Ideal supercell lattice (in Niggli form) of mapped configuration
    ///\param best_assignment[out]
    ///\parblock
    ///                   populated by the permutation of sites in the imported structure
    ///                   that maps them onto sites of the ideal crystal (excluding vacancies)
    ///\endparblock
    ///\param best_cost[in] optional parameter. Method will return false of no mapping is found better than 'best_cost'
    bool struc_to_configdof(const SimpleStructure &_struc,
                            MappedConfig &mapped_config,
                            Lattice &mapped_lat,
                            std::vector<Index> &best_assignment,
                            Eigen::Matrix3d &cart_op,
                            double best_cost = 1e20) const;


    ///\brief specify which lattices should be searched when mapping configurations
    void force_lattices(const std::vector<std::string> &lattice_names) const;

    ///\brief unset the enforcement of particular lattices (default behavior)
    void unforce_lattices() const;

    ///\brief specify to use restricted hermites when mapping
    void restricted() {
      m_restricted = true;
    }

    ///\brief returns true if lattices were set to be forced as candidates
    bool lattices_are_forced() const {
      return m_forced_superlat_map.size();
    }

  private:


    ///\brief Low-level routine to map a structure onto a ConfigDof if it is known to be ideal
    ///\param mapped_config[out] MappedConfig that is result of mapping procedure
    ///\param mapped_lat[out] Ideal supercell lattice (in Niggli form) of mapped configuration
    ///\param best_assignment[out]
    ///\parblock
    ///                   populated by the permutation of sites in the imported structure
    ///                   that maps them onto sites of the ideal crystal (excluding vacancies)
    ///\endparblock
    bool ideal_struc_to_configdof(const SimpleStructure &struc,
                                  MappedConfig &mapped_config,
                                  Lattice &mapped_lat,
                                  std::vector<Index> &best_assignment,
                                  Eigen::Matrix3d &cart_op) const;


    ///\brief Low-level routine to map a structure onto a ConfigDof. Does not assume structure is ideal
    ///\param mapped_config[out] MappedConfig that is result of mapping procedure
    ///\param mapped_lat[out] Ideal supercell lattice (in Niggli form) of mapped configuration
    ///\param best_assignment[out]
    ///\parblock
    ///                   populated by the permutation of sites in the imported structure
    ///                   that maps them onto sites of the ideal crystal (excluding vacancies)
    ///\endparblock
    ///\param best_cost[in] optional parameter. Method will return false of no mapping is found better than 'best_cost'
    bool deformed_struc_to_configdof(const SimpleStructure &_struc,
                                     MappedConfig &mapped_config,
                                     Lattice &mapped_lat,
                                     std::vector<Index> &best_assignment,
                                     Eigen::Matrix3d &cart_op,
                                     double best_cost = 1e20) const;

    ///\brief Low-level routine to map a structure onto a ConfigDof assuming a specific Lattice, without assuming structure is ideal
    ///       Will only identify mappings better than best_cost, and best_cost is updated to reflect cost of best mapping identified
    ///\param imposed_lat[in] Supercell Lattice onto which struc will be mapped
    ///\param best_cost Imposes an upper bound on cost of any mapping considered, and is updated to reflect best mapping encountered
    ///\param mapped_config[out] MappedConfig that is result of mapping procedure
    ///\param best_assignment[out]
    ///\parblock
    ///                   populated by the permutation of sites in the imported structure
    ///                   that maps them onto sites of the ideal crystal (excluding vacancies)
    ///\endparblock
    bool deformed_struc_to_configdof_of_lattice(const SimpleStructure &struc,
                                                const Lattice &imposed_lat,
                                                double &best_cost,
                                                MappedConfig &mapped_config,
                                                Lattice &mapped_lat,
                                                std::vector<Index> &best_assignment,
                                                Eigen::Matrix3d &cart_op) const;

  private:

    const PrimClex *m_pclex;
    ///Maps the supercell volume to a vector of Lattices with that volume
    mutable LatMapType m_superlat_map;
    mutable LatMapType m_forced_superlat_map;
    double m_lattice_weight;
    double m_max_volume_change;
    double m_min_va_frac;
    double m_max_va_frac;
    bool m_robust_flag, m_strict_flag, m_rotate_flag;
    bool m_restricted = false;
    double m_tol;
    std::vector<std::pair<std::string, Index> > m_fixed_components;
    const std::vector<Lattice> &_lattices_of_vol(Index prim_vol) const;
    const std::vector<Lattice> _lattices_of_vol_restricted(Index prim_vol, std::unordered_set<Eigen::Matrix3i, HermiteHash> &ref_hermites) const;
    bool hermite_adjacency(const Eigen::Matrix3i test, const Eigen::Matrix3i ref) const;
  };



  namespace ConfigMap_impl {
    // TODO: Refactor -- The two overloads of struc_to_configdof are *very* similar.
    // Assignment Problem Routines
    // Find cost matrix for displacements between ideal crystal and relaxed structure ('rstruc').
    // Returns false if 'rstruc' is incompatible with 'scel'
    bool calc_cost_matrix(const Supercell &scel,
                          const SimpleStructure &rstruc,
                          const Coordinate &trans,
                          const Eigen::Matrix3d &metric,
                          Eigen::MatrixXd &cost_matrix);

    // Assignment Problem Routines
    // Find cost matrix for displacements between ideal configuration and a relaxed structure ('rstruc').
    // Returns false if 'rstruc' is incompatible with 'scel'
    bool calc_cost_matrix(const Configuration &config,
                          const SimpleStructure &rstruc,
                          const Coordinate &trans,
                          const Eigen::Matrix3d &metric,
                          Eigen::MatrixXd &cost_matrix);



    // mapping routine. Return an ideal configuration corresponding to a relaxed structure.
    // Return false if 'rstruc' is incompatible with supercell (can happen frequently when vacancies are allowed)
    // Options:
    //   TRANSLATE = true -> rigid-translations are removed. (typically this option should be used, especially if you care about vacancies)
    //
    //   TRANSLATE = false -> rigid translations are not considered. (less robust but more efficient -- use only if you know rigid translations are small or zero)
    // TODO: Refactor -- The two overloads of struc_to_configdof are *very* similar.
    bool struc_to_configdof(const Supercell &scel,
                            SimpleStructure rstruc,
                            MappedConfig &mapped_config,
                            std::vector<Index> &best_assignment,
                            const bool translate_flag,
                            const double _tol);

    bool struc_to_configdof(const Configuration &config,
                            SimpleStructure rstruc,
                            MappedConfig &mapped_config,
                            std::vector<Index> &best_assignment,
                            const bool translate_flag,
                            const double _tol);

    // TODO: Refactor -- The two overloads of preconditioned_struc_to_configdof are *very* similar.
    /// Same as struc_to_configdof, except 'rstruc' is de-rotated and de-strained. Any deformation is instead specified by 'deformation'
    bool preconditioned_struc_to_configdof(const Supercell &scel,
                                           const SimpleStructure &rstruc,
                                           const Eigen::Matrix3d &deformation,
                                           MappedConfig &mapped_config,
                                           std::vector<Index> &best_assignment,
                                           const bool translate_flag,
                                           const double _tol);

    /// Same as struc_to_mapped_config, except 'rstruc' is de-rotated and de-strained. Any deformation is instead specified by 'deformation'
    bool preconditioned_struc_to_configdof(const Configuration &config,
                                           const SimpleStructure &rstruc,
                                           const Eigen::Matrix3d &deformation,
                                           MappedConfig &mapped_config,
                                           std::vector<Index> &best_assignment,
                                           const bool translate_flag,
                                           const double _tol);


  }

  namespace ConfigMapping {
    /// \brief Calculate the strain cost function of a MappedConfig using LatticeMap::calc_strain_cost()
    /// \param Nsites number of atoms in the relaxed structure, for proper normalization
    double strain_cost(double relaxed_lat_vol, const MappedConfig &mapped_config, Index Nsites);

    /// \brief Calculate the basis cost function of a MappedConfig as the mean-square displacement of its atoms
    /// \param Nsites number of atoms in the relaxed structure, for proper normalization
    double basis_cost(const MappedConfig &mapped_config, Index Nsites);


    /// \brief Find a supercell of prim_lat that is closest to relaxed_lat
    /// Considers all supercells having integer volumes on the range [min_vol, max_vol], wrt prim_lat.volume()
    /// Returns best lattice such that relaxed_lat = deformation * result * trans_mat
    /// \param sym_group point group of crystal (whose order is equal to or less than the lattice point group of prim_lat)
    /// used for supercell enumeration only
    /// \param trans_mat 3x3 integer matrix found by algorithm
    /// \param deformation 3x3 matrix of doubles found by algorithm
    /// \param best_cost Lattice mapping score, as defined in LatticeMap class
    Lattice find_nearest_super_lattice(const Lattice &prim_lat,
                                       const Lattice &relaxed_lat,
                                       const SymGroup &sym_group,
                                       Eigen::Matrix3d &trans_mat,
                                       Eigen::Matrix3d &deformation,
                                       double &best_cost,
                                       Index min_vol,
                                       Index max_vol,
                                       double _tol);

    /// \brief Find a supercell of prim_lat that is closest to relaxed_lat
    /// Considers all supercells contained in vector from_range
    /// Returns best lattice such that relaxed_lat = deformation * result * trans_mat
    /// \param sym_group point group of crystal (whose order is equal to or less than the lattice point group of prim_lat)
    /// used for supercell enumeration only
    /// \param trans_mat 3x3 integer matrix found by algorithm
    /// \param deformation 3x3 matrix of doubles found by algorithm
    /// \param best_cost Lattice mapping score, as defined in LatticeMap class
    Lattice find_nearest_super_lattice(const Lattice &prim_lat,
                                       const Lattice &relaxed_lat,
                                       const SymGroup &sym_group,
                                       Eigen::Matrix3d &trans_mat,
                                       Eigen::Matrix3d &deformation,
                                       double &best_cost,
                                       const std::vector<Lattice> &from_range,
                                       double _tol);


    ConfigMapperResult structure_mapping(Structure &host, Structure &other, double lattice_weight);
  }


}

#endif
