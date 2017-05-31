#ifndef ConfigDoF_HH
#define ConfigDoF_HH

#include <vector>
#include "casm/CASM_global_definitions.hh"
#include "casm/CASM_global_Eigen.hh"

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

    // typedefs to provide flexibility if we eventually change to a Eigen::Matrix<3,Eigen::Dynamic>
    typedef Eigen::MatrixXd displacement_matrix_t;

    // Can treat as a Eigen::VectorXd
    typedef displacement_matrix_t::ColXpr displacement_t;
    typedef displacement_matrix_t::ConstColXpr const_displacement_t;

    /// fixes alignment of m_deformation
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /// Initialize with number of sites -- defaults to zero
    ConfigDoF(Index N = 0, double _tol = TOL);

    ///Initialize with explicit occupation
    ConfigDoF(const std::vector<int> &_occupation, double _tol = TOL);

    ///Initialize with explicit properties
    ConfigDoF(const std::vector<int> &_occupation, const Eigen::MatrixXd &_displacement, const Eigen::Matrix3d &_deformation, double _tol = TOL);


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


    // -- Displacement ------------------

    displacement_t disp(Index i) {
      return m_displacement.col(i);
    };

    const_displacement_t disp(Index i) const {
      return m_displacement.col(i);
    };

    /// set_displacement ensures that ConfigDoF::size() is compatible with
    /// _displacement.cols() (i.e., number of sites)
    /// or if ConfigDoF::size()==0, sets ConfigDoF::size() to _displacement.cols()
    void set_displacement(const displacement_matrix_t &_displacement);

    displacement_matrix_t &displacement() {
      return m_displacement;
    };

    const displacement_matrix_t &displacement() const {
      return m_displacement;
    };

    bool has_displacement() const {
      return size() != 0 && displacement().cols() == size();
    }

    void clear_displacement();


    // -- Deformation ------------------

    /// set_deformation sets ConfigDoF::has_deformation() to true
    void set_deformation(const Eigen::Matrix3d &_deformation);

    Eigen::Matrix3d &deformation() {
      return m_deformation;
    }

    const Eigen::Matrix3d &deformation() const {
      return m_deformation;
    }

    const double &F(Index i, Index j) const {
      return m_deformation(i, j);
    }

    double &F(Index i, Index j) {
      assert(has_deformation() && "Non-const method ConfigDoF::F() should only be called after ConfigDoF::set_deformation()!!!");
      return m_deformation(i, j);
    }

    bool has_deformation() const {
      return m_has_deformation;
    }

    void clear_deformation();


    // -- Specie ID ------------------

    /// Access specie id data
    std::vector<std::vector<Index> > &specie_id() {
      return m_specie_id;
    }

    const std::vector<std::vector<Index> > &specie_id() const {
      return m_specie_id;
    }

    bool has_specie_id() const {
      return size() != 0 && specie_id().size() == size();
    }

    void clear_specie_id();


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
    ///   prim().basis[ sublat(i) ].site_occupant[ occupation[i]] -> Molecule on site i
    ///   This means that for the background structure, 'occupation' is all 0
    ///
    /// Configuration sites are arranged by basis, and then prim:
    ///   occupation: [basis0                |basis1               |basis2          |...] up to prim.basis.size()
    ///       basis0: [prim0|prim1|prim2|...] up to supercell.volume()
    ///
    ///
    std::vector<int> m_occupation;

    /// A VectorXd for each site in the Configuration to describe displacements condensed in matrix form  -- This a 3xN matrix whose columns are the displacement of
    /// each of the N sites of the configuration
    displacement_matrix_t m_displacement;

    ///Describes possible strains that may have been applied to the Configuration -- This is the matrix that relates the reference lattice vectors
    ///to the deformed lattice vectors via L_deformed = m_deformation * L_reference -- (L is a 3x3 matrix whose columns are the lattice vectors)
    Eigen::Matrix3d m_deformation;

    bool m_has_deformation;

    /// Use to track specie ID during KMC, vector at each site to handle atoms
    /// in a molecule
    std::vector<std::vector<Index> > m_specie_id;

    /// Tolerance used for transformation to canonical form -- used also for comparisons, since
    /// Since comparisons are only meaningful to within the tolerance used for finding the canonical form
    /// (This is relevant only for displacement and deformation degrees of freedom
    mutable double m_tol;

  };

  jsonParser &to_json(const ConfigDoF &value, jsonParser &json);

  void from_json(ConfigDoF &value, const jsonParser &json);

  // Calculate transformed ConfigDoF from PermuteIterator via
  //   apply(permute_iterator, dof)
  ConfigDoF &apply(const PermuteIterator &it, ConfigDoF &dof);

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
