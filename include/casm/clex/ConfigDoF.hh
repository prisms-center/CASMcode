#ifndef ConfigDoF_HH
#define ConfigDoF_HH

#include "casm/container/Array.hh"
//#include "casm/symmetry/PermuteIterator.hh"

namespace CASM {


  class PermuteIterator;
  typedef Array<double> Correlation;
  class PrimClex;
  class Supercell;
  class Clexulator;


  /**
   * This is just a container class for the different degrees of freedom a Configuration
   * might have. Contains an id, an Array<int> that tells you the current occupant of each
   * site, an Eigen::MatrixXd that tells you the displacements at each site, and a LatticeStrain
   * that tells you the strain of the Configuration. Everything is public.
   */

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
    ConfigDoF(Index N = 0);

    ///Initialize with explicit occupation
    ConfigDoF(const Array<int> &_occupation);

    ///Initialize with explicit properties
    ConfigDoF(const Array<int> &_occupation, const Eigen::MatrixXd &_displacement, const Eigen::Matrix3d &_deformation);

  public:

    /// *** ACCESSORS ***
    Index size() const {
      return m_N;
    }

    bool operator==(const ConfigDoF &RHS) const;

    int &occ(Index i) {
      return m_occupation[i];
    };

    const int &occ(Index i) const {
      return m_occupation[i];
    };

    const Array<int> &occupation() const {
      return m_occupation;
    };

    displacement_t disp(Index i) {
      return m_displacement.col(i);
    };

    const_displacement_t disp(Index i) const {
      return m_displacement.col(i);
    };

    const displacement_matrix_t &displacement() const {
      return m_displacement;
    };

    const Eigen::Matrix3d &deformation() const {
      return m_deformation;
    }

    const double &F(Index i, Index j) const {
      return m_deformation(i, j);
    }

    double &F(Index i, Index j) {
      assert(is_strained() && "Non-const method ConfigDoF::F() should only be called after ConfigDoF::set_deformation()!!!");
      return m_deformation(i, j);
    }

    bool is_strained() const {
      return m_is_strained;
    }

    bool has_displacement() const {
      return size() != 0 && displacement().cols() == size();
    }

    bool has_occupation() const {
      return size() != 0 && occupation().size() == size();
    }

    //**** MUTATORS ****
    /// set_occupation ensures that ConfigDoF::size() is compatible with _occupation.size()
    /// or if ConfigDoF::size()==0, sets ConfigDoF::size() to _occupation.size()
    void set_occupation(const Array<int> &_occupation);

    /// set_displacement ensures that ConfigDoF::size() is compatible with _displacement.cols() (i.e., number of sites)
    /// or if ConfigDoF::size()==0, sets ConfigDoF::size() to _displacement.cols()
    void set_displacement(const displacement_matrix_t &_displacement);

    /// set_deformation sets ConfigDoF::is_strained() to true
    void set_deformation(const Eigen::Matrix3d &_deformation);

    void clear();

    void swap(ConfigDoF &RHS);

    // pass iterator, 'it_begin', to first translation operation (indexed 0,0) to determine if configdof is a primitive cell
    bool is_primitive(PermuteIterator it_begin, double tol = TOL) const;

    ReturnArray<PermuteIterator> factor_group(PermuteIterator it_begin, PermuteIterator it_end, double tol = TOL) const;

    bool is_canonical(PermuteIterator it_begin, PermuteIterator it_end, double tol = TOL) const;

    bool is_canonical(PermuteIterator it_begin, PermuteIterator it_end, Array<PermuteIterator> &factor_group, double tol = TOL) const;

    ConfigDoF canonical_form(PermuteIterator it_begin, PermuteIterator it_end, double tol = TOL) const;

    ConfigDoF canonical_form(PermuteIterator it_begin, PermuteIterator it_end, PermuteIterator &it_canon, double tol = TOL) const;

    ConfigDoF canonical_form(PermuteIterator it_begin, PermuteIterator it_end, PermuteIterator &it_canon, Array<PermuteIterator> &factor_group, double tol = TOL) const;

    //**** I/O ****
    jsonParser &to_json(jsonParser &json) const;
    void from_json(const jsonParser &json);


  private:

    // ***DON'T FORGET: If you add something here, also update ConfigDoF::swap!!

    
    /// \brief Number of sites in the Configuration
    Index m_N;

    ///With one value for each site in the Configuration, this Array describes which occupant is at each of the 'N' sites of the configuration
    Array<int> m_occupation;

    /// A VectorXd for each site in the Configuration to describe displacements condensed in matrix form  -- This a 3xN matrix whose columns are the displacement of
    /// each of the N sites of the configuration
    displacement_matrix_t m_displacement;

    ///Describes possible strains that may have been applied to the Configuration -- This is the matrix that relates the reference lattice vectors
    ///to the deformed lattice vectors via L_deformed = m_deformation * L_reference -- (L is a 3x3 matrix whose columns are the lattice vectors)
    Eigen::Matrix3d m_deformation;

    bool m_is_strained;

    /// Tolerance used for transformation to canonical form -- used also for comparisons, since
    /// Since comparisons are only meaningful to within the tolerance used for finding the canonical form
    /// (This is relevant only for displacement and deformation degrees of freedom
    mutable double m_tol;

    // *** PRIVATE NON-CONST ACCESS
    Array<int> &_occupation() {
      return m_occupation;
    };

    displacement_matrix_t &_displacement() {
      return m_displacement;
    };

    Eigen::Matrix3d &_deformation() {
      return m_deformation;
    }

    // *** private implementation of is_canonical() and canonical_form()
    bool _is_canonical(PermuteIterator it_begin, PermuteIterator it_end, Array<PermuteIterator> *fg_ptr = NULL, double tol = TOL) const;

    ConfigDoF _canonical_form(PermuteIterator it_begin, PermuteIterator it_end, PermuteIterator &it_canon, Array<PermuteIterator> *fg_ptr = NULL, double tol = TOL) const;

  };

  jsonParser &to_json(const ConfigDoF &value, jsonParser &json);

  void from_json(ConfigDoF &value, const jsonParser &json);


  // Calculate transformed ConfigDoF from PermuteIterator via
  //   transformed_configdof = permute_iterator * configdof
  ConfigDoF operator*(const PermuteIterator &it, const ConfigDoF &dof);

  void swap(ConfigDoF &A, ConfigDoF &B);

  /// \brief Returns correlations using 'clexulator'. Supercell needs a correctly populated neighbor list.
  Correlation correlations(const ConfigDoF &configdof, const Supercell &scel, Clexulator &clexulator);
  
  /// \brief Returns correlations using 'clexulator'. Supercell needs a correctly populated neighbor list.
  Eigen::VectorXd correlations_vec(const ConfigDoF &configdof, const Supercell &scel, Clexulator &clexulator);
  
  /// \brief Returns num_each_molecule[ molecule_type], where 'molecule_type' is ordered as Structure::get_struc_molecule()
  ReturnArray<int> get_num_each_molecule(const ConfigDoF &configdof, const Supercell &scel); 
  
  /// \brief Returns num_each_molecule(molecule_type), where 'molecule_type' is ordered as Structure::get_struc_molecule()
  Eigen::VectorXi get_num_each_molecule_vec(const ConfigDoF &configdof, const Supercell &scel);
  
  /// \brief Returns comp_n, the number of each molecule per primitive cell, ordered as Structure::get_struc_molecule()
  Eigen::VectorXd comp_n(const ConfigDoF &configdof, const Supercell &scel);
  
  /// \brief Return a super ConfigDoF (occupation only)
  ConfigDoF super_configdof_occ(PrimClex &primclex, const Eigen::Matrix3i& transf_mat, const std::string &motif_configname);
  

}

#endif
