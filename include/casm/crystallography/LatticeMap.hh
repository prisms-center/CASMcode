#ifndef LATTICEMAP_HH
#define LATTICEMAP_HH

#include "casm/container/Counter.hh"
#include "casm/crystallography/Lattice.hh"

namespace CASM {

  /** \ingroup Lattice
   *  @{
   */

  /// Find the ideal mapping of Lattice _ideal onto Lattice _strained
  /// Denoting _ideal.lat_column_mat() as 'L1' and _strained.lat_column_mat() as 'L2', we want a mapping
  ///            L2 = F*L1*N
  /// where F is an arbitrary 3x3 matrix and 'N' is a 3x3 unimodular (i.e., determinant=+/-1) integer matrix
  /// such that cost(F) is minimized with respect to the matrix 'N'
  /// For the cost function we use:
  ///         cost(F) = (_strained.volume()/num_atoms)^(2/3) * trace(D.transpose()*D) / g
  /// where 'D' is the isovolumetric strain
  ///         D = F/det(F)^(1/3)-Identity
  /// and 'g' is a const geometric factor.  The cost function corresponds to the mean-square displacement of a point
  /// on the surface of a sphere having V=_strained.volume()/num_atoms (i.e., the atomic volume of the strained crystal)
  /// when the sphere is deformed at constant volume by F/det(F)^(1/3)


  class LatticeMap {
  public:
    typedef Eigen::Matrix<double, 3, 3, Eigen::DontAlign> DMatType;
    typedef Eigen::Matrix<int, 3, 3, Eigen::DontAlign> IMatType;

    LatticeMap(Lattice const &_ideal,
               Lattice const &_strained,
               Index _num_atoms,
               double _tol /*= TOL*/,
               int _range /*= 2*/,
               std::vector<SymOp> const &_point_group/*={}*/);

    LatticeMap(Eigen::Ref<const DMatType> const &_ideal,
               Eigen::Ref<const DMatType> const &_strained,
               Index _num_atoms,
               double _tol /*= TOL*/,
               int _range /*= 2*/,
               std::vector<SymOp> const &_point_group/*={}*/);

    void reset();

    // Finds the smallest strain tensor (in terms of Frobenius norm) that deforms (*this) into a lattice symmetrically equivalent to 'strained_lattice'
    LatticeMap const &best_strain_mapping() const;

    LatticeMap const &next_mapping_better_than(double max_cost) const;

    double strain_cost() const {
      return m_cost;
    }

    const DMatType &matrixN() const {
      return m_N;
    }

    const DMatType &matrixF() const {
      return m_F;
    }

    // calculated strain cost function given deformation gradient 'F' and volume of the relaxed cell
    static double calc_strain_cost(const Eigen::Matrix3d &F, double relaxed_atomic_vol);

  private:
    DMatType m_L1, m_L2;
    //Conversion matrices:
    //  m_N = m_U * m_inv_count().inverse() * m_V_inv
    DMatType m_U, m_V_inv;

    // m_scale = (det(m_L2)/det(m_L1))^(2/3) = det(m_F)^(2/3)
    double m_scale;
    // m_atomic_factor = (det(m_L2)/num_atoms)^(2/3)
    double m_atomic_factor;
    double m_tol;
    int m_range;
    std::vector<Eigen::Matrix3i> const *m_mvec_ptr;
    std::vector<Eigen::Matrix3i> m_fsym_mats;
    mutable double m_cost;
    mutable Index m_currmat;
    mutable DMatType m_F, m_N, m_cache;


    Eigen::Matrix3i const &inv_mat() const {
      return (*m_mvec_ptr)[m_currmat];
    }

    Index n_mat() const {
      return m_mvec_ptr->size();
    }

    bool _check_canonical() const;

    LatticeMap const &_next_mapping_better_than(double max_cost) const;

    // use m_F and m_atomic_vol to calculate strain cost
    double _calc_strain_cost() const;

  };

  /** @} */
}
#endif
