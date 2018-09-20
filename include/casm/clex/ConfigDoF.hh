#ifndef ConfigDoF_HH
#define ConfigDoF_HH

#include <vector>
#include "casm/CASM_global_definitions.hh"
#include "casm/CASM_global_Eigen.hh"
#include "casm/clex/ConfigDoFValues.hh"
namespace CASM {


  class PermuteIterator;
  class PrimClex;
  class Supercell;
  class Clexulator;
  class jsonParser;

  /// \brief A container class for the different degrees of freedom a Configuration
  /// might have
  ///
  /// Contains an std::vector<int> that tells you the current occupant of each
  /// site, an Eigen::MatrixXd that tells you the displacements at each site, and
  /// a LatticeStrain that tells you the strain of the Configuration. Everything
  /// is public.
  ///
  /// \ingroup Configuration
  ///
  class ConfigDoF {

  public:

    // Can treat as a Eigen::VectorXd
    //typedef displacement_matrix_t::ColXpr displacement_t;
    //typedef displacement_matrix_t::ConstColXpr const_displacement_t;

    /// Initialize with number of sites -- defaults to zero
    ConfigDoF(Index N = 0, double _tol = TOL);

    ///Initialize with explicit occupation
    ConfigDoF(const std::vector<int> &_occupation, double _tol = TOL);

    Index size() const {
      return m_N;
    }

    double tol() const {
      return m_tol;
    }

    void clear();


    // -- Occupation ------------------

    int &occ(Index i) {
      return m_occupation[i];
    };

    const int &occ(Index i) const {
      return m_occupation[i];
    };

    /// set_occupation ensures that ConfigDoF::size() is compatible with _occupation.size()
    /// or if ConfigDoF::size()==0, sets ConfigDoF::size() to _occupation.size()
    void set_occupation(const std::vector<int> &_occupation);

    std::vector<int> &occupation() {
      return m_occupation;
    };

    const std::vector<int> &occupation() const {
      return m_occupation;
    };

    bool has_occupation() const {
      return size() != 0 && occupation().size() == size();
    }

    void clear_occupation();


    bool has_global_dof(DoFKey const &_key);
    void reset_global_dof(DoFKey const &_key);
    void set_global_dof(DoFKey const &_key, Eigen::Ref<const Eigen::VectorXd> const &_val);

    bool has_local_dof(DoFKey const &_key);
    void reset_local_dof(DoFKey const &_key);
    void set_local_dof(DoFKey const &_key, Eigen::Ref<const Eigen::MatrixXd> const &_val);

    ConfigDoF &apply_sym(const PermuteIterator &it);

    void swap(ConfigDoF &RHS);

    //**** I/O ****
    jsonParser &to_json(jsonParser &json) const;
    void from_json(const jsonParser &json);


  private:

    // ***DON'T FORGET: If you add something here, also update ConfigDoF::swap!!


    /// \brief Number of sites in the Configuration
    Index m_N;

    /// With one value for each site in the Configuration, this std::vector
    /// describes which occupant is at each of the 'N' sites of the configuration
    ///
    /// 'occupation' is a list of the indices describing the occupants in each crystal site.
    ///   prim().basis()[ sublat(i) ].site_occupant[ occupation[i]] -> Molecule on site i
    ///   This means that for the background structure, 'occupation' is all 0
    ///
    /// Configuration sites are arranged by basis, and then prim:
    ///   occupation: [basis0                |basis1               |basis2          |...] up to prim.basis().size()
    ///       basis0: [prim0|prim1|prim2|...] up to supercell.volume()
    ///
    ///
    std::vector<int> m_occupation;

    std::vector<GlobalContinuousConfigDoFValues> m_global_vals;

    std::vector<LocalContinuousConfigDoFValues> m_local_vals;

    /// Tolerance used for transformation to canonical form -- used also for comparisons, since
    /// Since comparisons are only meaningful to within the tolerance used for finding the canonical form
    /// (This is relevant only for displacement and deformation degrees of freedom
    mutable double m_tol;

  };

  jsonParser &to_json(const ConfigDoF &value, jsonParser &json);

  void from_json(ConfigDoF &value, const jsonParser &json);

  void swap(ConfigDoF &A, ConfigDoF &B);

  /// \brief Returns correlations using 'clexulator'. Supercell needs a correctly populated neighbor list.
  Eigen::VectorXd correlations(const ConfigDoF &configdof, const Supercell &scel, Clexulator &clexulator);

  /// \brief Returns num_each_molecule(molecule_type), where 'molecule_type' is ordered as Structure::get_struc_molecule()
  Eigen::VectorXi num_each_molecule(const ConfigDoF &configdof, const Supercell &scel);

  /// \brief Returns comp_n, the number of each molecule per primitive cell, ordered as Structure::get_struc_molecule()
  Eigen::VectorXd comp_n(const ConfigDoF &configdof, const Supercell &scel);

  inline
  void reset_properties(ConfigDoF &_dof) {
    return;
  }

}

#endif
