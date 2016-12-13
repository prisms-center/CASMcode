#ifndef CONFIGMAPPING_HH
#define CONFIGMAPPING_HH
#include "casm/CASM_global_definitions.hh"
#include <vector>
namespace CASM {
  class Supercell;
  class Lattice;
  class SymGroup;
  class PrimClex;
  class Configuration;
  class ConfigDoF;

  /// A class for mapping an arbitrary crystal structure as a configuration of a crystal template
  /// as described by a PrimClex.  ConfigMapper manages options for the mapping algorithm and mapping cost function
  /// It also caches some information about supercell lattices so that batch imports are more efficient
  ///
  /// \ingroup Configuration
  class ConfigMapper {
  public:
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
    ConfigMapper(PrimClex &_pclex,
                 double _lattice_weight,
                 double _max_volume_change = 0.5,
                 int _options = robust, // this should actually be a bitwise-OR of ConfigMapper::Options
                 double _tol = TOL);


    PrimClex &primclex() const {
      return *m_pclex;
    }

    void set_primclex(PrimClex &_pclex) {
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
    ///\param imported_name is populated by the configname given to the imported structure
    ///                     (or an existing equivalent structure if one exists)
    ///\param relaxation_properties is populated by relaxation properties:
    ///                     'lattice_deformation':    lattice mapping score
    ///                     'basis_deformation':      atomic mapping score
    ///                     'volume_relaxation':      V/V_ideal
    ///                     'relaxation_deformation': 3x3 tensor describing cell relaxation
    ///\param best_assignment is populated by the permutation of sites in the imported structure
    ///                       that maps them onto sites of the ideal crystal (excluding vacancies)
    ///\param cart_op is populated by the cartesian isometry that rotates the imported structure
    ///               onto the coordinate system of the ideal crystal
    bool import_structure_occupation(const fs::path &pos_path,
                                     std::string &imported_name,
                                     jsonParser &relaxation_properties,
                                     std::vector<Index> &best_assignment,
                                     Eigen::Matrix3d &cart_op) const;

    ///\brief imports structure specified by '_struc' into primclex()
    /// \param update_struc if true, applies lattice similarity and/or rotation to _struc. If false, _struc is unchanged
    /// update_struc=true is useful for preserving mapping info that is lost when only the occupation is imported
    bool import_structure_occupation(BasicStructure<Site> &_struc,
                                     std::string &imported_name,
                                     jsonParser &relaxation_properties,
                                     std::vector<Index> &best_assignment,
                                     Eigen::Matrix3d &cart_op,
                                     bool update_struc = false) const;

    ///\brief imports structure specified by '_struc' into primclex()
    ///\param hint_ptr[in]
    ///\parblock
    ///                provides a suggestion for which Configuration _struc should map onto
    ///                The hint is used to reduce search times, but may be used in the future
    ///                in combination with Option 'strict' to force mapping onto a particular configuration
    ///                or be used to provide user reports of the form "Suggested mapping: 0.372; Optimal mapping: 0.002"
    ///\endparblock
    /// \param update_struc if true, applies lattice similarity and/or rotation to _struc. If false, _struc is unchanged
    /// update_struc=true is useful for preserving mapping info that is lost when only the occupation is imported
    bool import_structure_occupation(BasicStructure<Site> &_struc,
                                     const Configuration *hint_ptr,
                                     std::string &imported_name,
                                     jsonParser &relaxation_properties,
                                     std::vector<Index> &best_assignment,
                                     Eigen::Matrix3d &cart_op,
                                     bool update_struc = false) const;


    ///\brief imports structure specified by 'pos_path' into primclex() by finding optimal mapping
    ///       unlike import_structure_occupation, displacements and strain are preserved
    ///
    ///\param imported_name[out]
    ///\parblock
    ///                 populated by the configname given to the imported structure
    ///                 (or an existing equivalent structure if one exists)
    ///\endparblock
    ///
    ///\param relaxation_properties[out]
    ///\parblock
    ///                 populated by relaxation properties:
    ///                     'lattice_deformation':    lattice mapping score
    ///                     'basis_deformation':      atomic mapping score
    ///                     'volume_change':          V/V_ideal
    ///                     'relaxation_deformation': Not recorded in this case, because it is a DoF
    ///\endparblock
    ///
    ///\param best_assignment[out]
    ///\parblock
    ///                   populated by the permutation of sites in the imported structure
    ///                   that maps them onto sites of the ideal crystal (excluding vacancies)
    ///\endparblock
    ///
    ///\param cart_op[out]
    ///\parblock
    ///                   populated by the cartesian isometry that rotates the imported structure
    ///                   onto the coordinate system of the ideal crystal
    ///\endparblock
    bool import_structure(const fs::path &pos_path,
                          std::string &imported_name,
                          jsonParser &relaxation_properties,
                          std::vector<Index> &best_assignment,
                          Eigen::Matrix3d &cart_op) const;

    ///\brief imports structure specified by '_struc' into primclex() by finding optimal mapping
    ///       unlike import_structure_occupation, displacements and strain are preserved
    bool import_structure(const BasicStructure<Site> &_struc,
                          std::string &imported_name,
                          jsonParser &relaxation_properties,
                          std::vector<Index> &best_assignment,
                          Eigen::Matrix3d &cart_op) const;

    ///\brief Low-level routine to map a structure onto a ConfigDof
    ///\param mapped_configdof[out] ConfigDoF that is result of mapping procedure
    ///\param mapped_lat[out] Ideal supercell lattice (in Niggli form) of mapped configuration
    bool struc_to_configdof(const BasicStructure<Site> &_struc,
                            ConfigDoF &mapped_configdof,
                            Lattice &mapped_lat) const;

    ///\brief Low-level routine to map a structure onto a ConfigDof
    ///\param mapped_configdof[out] ConfigDoF that is result of mapping procedure
    ///\param mapped_lat[out] Ideal supercell lattice (in Niggli form) of mapped configuration
    ///\param best_assignment[out]
    ///\parblock
    ///                   populated by the permutation of sites in the imported structure
    ///                   that maps them onto sites of the ideal crystal (excluding vacancies)
    ///\endparblock
    ///\param best_cost[in] optional parameter. Method will return false of no mapping is found better than 'best_cost'
    bool struc_to_configdof(const BasicStructure<Site> &_struc,
                            ConfigDoF &mapped_configdof,
                            Lattice &mapped_lat,
                            std::vector<Index> &best_assignment,
                            Eigen::Matrix3d &cart_op,
                            double best_cost = 1e20) const;


    ///\brief Low-level routine to map a structure onto a ConfigDof if it is known to be ideal
    ///\param mapped_configdof[out] ConfigDoF that is result of mapping procedure
    ///\param mapped_lat[out] Ideal supercell lattice (in Niggli form) of mapped configuration
    ///\param best_assignment[out]
    ///\parblock
    ///                   populated by the permutation of sites in the imported structure
    ///                   that maps them onto sites of the ideal crystal (excluding vacancies)
    ///\endparblock
    bool ideal_struc_to_configdof(const BasicStructure<Site> &struc,
                                  ConfigDoF &mapped_config_dof,
                                  Lattice &mapped_lat,
                                  std::vector<Index> &best_assignment,
                                  Eigen::Matrix3d &cart_op) const;


    ///\brief Low-level routine to map a structure onto a ConfigDof. Does not assume structure is ideal
    ///\param mapped_configdof[out] ConfigDoF that is result of mapping procedure
    ///\param mapped_lat[out] Ideal supercell lattice (in Niggli form) of mapped configuration
    ///\param best_assignment[out]
    ///\parblock
    ///                   populated by the permutation of sites in the imported structure
    ///                   that maps them onto sites of the ideal crystal (excluding vacancies)
    ///\endparblock
    ///\param best_cost[in] optional parameter. Method will return false of no mapping is found better than 'best_cost'
    bool deformed_struc_to_configdof(const BasicStructure<Site> &_struc,
                                     ConfigDoF &mapped_config_dof,
                                     Lattice &mapped_lat,
                                     std::vector<Index> &best_assignment,
                                     Eigen::Matrix3d &cart_op,
                                     double best_cost = 1e20) const;

    ///\brief Low-level routine to map a structure onto a ConfigDof assuming a specific Lattice, without assuming structure is ideal
    ///       Will only identify mappings better than best_cost, and best_cost is updated to reflect cost of best mapping identified
    ///\param imposed_lat[in] Supercell Lattice onto which struc will be mapped
    ///\param best_cost Imposes an upper bound on cost of any mapping considered, and is updated to reflect best mapping encountered
    ///\param mapped_configdof[out] ConfigDoF that is result of mapping procedure
    ///\param best_assignment[out]
    ///\parblock
    ///                   populated by the permutation of sites in the imported structure
    ///                   that maps them onto sites of the ideal crystal (excluding vacancies)
    ///\endparblock
    bool deformed_struc_to_configdof_of_lattice(const BasicStructure<Site> &struc,
                                                const Lattice &imposed_lat,
                                                double &best_cost,
                                                ConfigDoF &mapped_configdof,
                                                Lattice &mapped_lat,
                                                std::vector<Index> &best_assignment,
                                                Eigen::Matrix3d &cart_op) const;

  private:
    PrimClex *m_pclex;
    mutable std::map<Index, std::vector<Lattice> > m_superlat_map;
    double m_lattice_weight;
    double m_max_volume_change;
    double m_min_va_frac;
    double m_max_va_frac;
    bool m_robust_flag, m_strict_flag, m_rotate_flag;
    double m_tol;
    std::vector<std::pair<std::string, Index> > m_fixed_components;
    const std::vector<Lattice> &_lattices_of_vol(Index prim_vol) const;
  };

  namespace ConfigMap_impl {

    // Assignment Problem Routines
    // Find cost matrix for displacements between ideal crystal and relaxed structure ('rstruc').
    // Returns false if 'rstruc' is incompatible with 'scel'
    bool calc_cost_matrix(const Supercell &scel,
                          const BasicStructure<Site> &rstruc,
                          const Coordinate &trans,
                          const Eigen::Matrix3d &metric,
                          Eigen::MatrixXd &cost_matrix);

    // Assignment Problem Routines
    // Find cost matrix for displacements between ideal configuration and a relaxed structure ('rstruc').
    // Returns false if 'rstruc' is incompatible with 'scel'
    bool calc_cost_matrix(const Configuration &config,
                          const BasicStructure<Site> &rstruc,
                          const Coordinate &trans,
                          const Eigen::Matrix3d &metric,
                          Eigen::MatrixXd &cost_matrix);

    //\JSB



    // mapping routine. Return an ideal configuration corresponding to a relaxed structure.
    // Return false if 'rstruc' is incompatible with supercell (can happen frequently when vacancies are allowed)
    // Options:
    //   TRANSLATE = true -> rigid-translations are removed. (typically this option should be used, especially if you care about vacancies)
    //
    //   TRANSLATE = false -> rigid translations are not considered. (less robust but more efficient -- use only if you know rigid translations are small or zero)

    bool struc_to_configdof(const Supercell &scel,
                            BasicStructure<Site> rstruc,
                            ConfigDoF &config_dof,
                            std::vector<Index> &best_assignment,
                            const bool translate_flag,
                            const double _tol);

    /// Same as struc_to_configdof, except 'rstruc' is de-rotated and de-strained. Any deformation is instead specified by 'deformation'
    bool preconditioned_struc_to_configdof(const Supercell &scel,
                                           const BasicStructure<Site> &rstruc,
                                           const Eigen::Matrix3d &deformation,
                                           ConfigDoF &config_dof,
                                           std::vector<Index> &best_assignment,
                                           const bool translate_flag,
                                           const double _tol);

    bool struc_to_configdof(const Configuration &config,
                            BasicStructure<Site> rstruc,
                            ConfigDoF &config_dof,
                            std::vector<Index> &best_assignment,
                            const bool translate_flag,
                            const double _tol);

    /// Same as struc_to_configdof, except 'rstruc' is de-rotated and de-strained. Any deformation is instead specified by 'deformation'
    bool preconditioned_struc_to_configdof(const Configuration &config,
                                           const BasicStructure<Site> &rstruc,
                                           const Eigen::Matrix3d &deformation,
                                           ConfigDoF &config_dof,
                                           std::vector<Index> &best_assignment,
                                           const bool translate_flag,
                                           const double _tol);


  }

  namespace ConfigMapping {
    /// \brief Calculate the strain cost function of a ConfigDoF using LatticeMap::calc_strain_cost()
    /// \param Nsites number of atoms in the relaxed structure, for proper normalization
    double strain_cost(const Lattice &relaxed_lat, const ConfigDoF &_dof, Index Nsites);

    /// \brief Calculate the basis cost function of a ConfigDoF as the mean-square displacement of its atoms
    /// \param Nsites number of atoms in the relaxed structure, for proper normalization
    double basis_cost(const ConfigDoF &_dof, Index Nsites);


    Lattice find_nearest_super_lattice(const Lattice &prim_lat,
                                       const Lattice &relaxed_lat,
                                       const SymGroup &sym_group,
                                       Eigen::Matrix3d &trans_mat,
                                       Eigen::Matrix3d &deformation,
                                       double &best_cost,
                                       Index min_vol,
                                       Index max_vol,
                                       double _tol);

    Lattice find_nearest_super_lattice(const Lattice &prim_lat,
                                       const Lattice &relaxed_lat,
                                       const SymGroup &sym_group,
                                       Eigen::Matrix3d &trans_mat,
                                       Eigen::Matrix3d &deformation,
                                       double &best_cost,
                                       const std::vector<Lattice> &from_range,
                                       double _tol);


  }


}

#endif
