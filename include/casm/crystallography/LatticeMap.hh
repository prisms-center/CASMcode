#ifndef LATTICEMAP_HH
#define LATTICEMAP_HH

#include "casm/container/Counter.hh"
#include "casm/crystallography/Lattice.hh"

namespace CASM {
  namespace xtal {

    /** \ingroup Lattice
     *  @{
     */


    class StrainCostCalculator {
    public:
      StrainCostCalculator(Eigen::Ref<const Eigen::MatrixXd> const &strain_gram_mat = Eigen::MatrixXd::Identity(9, 9));

      //\brief Isotropic strain cost, without gram matrix
      static double iso_strain_cost(Eigen::Matrix3d const &F, double relaxed_atomic_vol);

      //\brief Isotropic strain cost, without gram matrix
      static double iso_strain_cost(Eigen::Matrix3d const &F, double relaxed_atomic_vol, double _vol_factor);

      // \brief Volumetric factor : pow(abs(F.determinant()),1./3.), used to normalize the strain cost to make it volume-independent
      static double vol_factor(Eigen::Matrix3d const &F) {
        return pow(std::abs(F.determinant()), 1. / 3.);
      }

      //\brief Anisotropic strain cost; utilizes stored gram matrix to compute strain cost
      double strain_cost(Eigen::Matrix3d const &F, double relaxed_atomic_vol)const;

      //\brief Anisotropic strain cost; utilizes stored gram matrix to compute strain cost
      double strain_cost(Eigen::Matrix3d const &F, double relaxed_atomic_vol, double _vol_factor)const;


    private:
      Eigen::MatrixXd m_gram_mat;
      bool m_sym_cost;

      mutable Eigen::Matrix3d m_cache;
      mutable Eigen::Matrix3d m_cache_inv;
    };

    /// Find the parent mapping of Lattice _parent onto Lattice _child
    /// Denoting _parent.lat_column_mat() as 'parent' and _child.lat_column_mat() as 'child', we want a mapping
    ///            child = F*parent*N
    /// where F is an arbitrary 3x3 matrix and 'N' is a 3x3 unimodular (i.e., determinant=+/-1) integer matrix
    /// such that cost(F) is minimized with respect to the matrix 'N'
    /// For the cost function we use:
    ///         cost(F) = (_child.volume()/num_atoms)^(2/3) * trace(D.transpose()*D) / g
    /// where 'D' is the isovolumetric strain
    ///         D = F/det(F)^(1/3)-Identity
    /// and 'g' is a const geometric factor.  The cost function corresponds to the mean-square displacement of a point
    /// on the surface of a sphere having V=_child.volume()/num_atoms (i.e., the atomic volume of the child crystal)
    /// when the sphere is deformed at constant volume by F/det(F)^(1/3)


    class LatticeMap {
    public:
      typedef Eigen::Matrix<double, 3, 3, Eigen::DontAlign> DMatType;
      typedef Eigen::Matrix<int, 3, 3, Eigen::DontAlign> IMatType;

      LatticeMap(Lattice const &_parent,
                 Lattice const &_child,
                 Index _num_atoms,
                 double _tol,
                 int _range /*= 2*/,
                 std::vector<CASM::SymOp> const &_point_group/*={}*/,
                 Eigen::Ref<const Eigen::MatrixXd> const &strain_gram_mat = Eigen::MatrixXd::Identity(9, 9),
                 double _init_better_than = 1e20);

      LatticeMap(Eigen::Ref<const DMatType> const &_parent,
                 Eigen::Ref<const DMatType> const &_child,
                 Index _num_atoms,
                 double _tol,
                 int _range /*= 2*/,
                 std::vector<CASM::SymOp> const &_point_group/*={}*/,
                 Eigen::Ref<const Eigen::MatrixXd> const &strain_gram_mat = Eigen::MatrixXd::Identity(9, 9),
                 double _init_better_than = 1e20);

      void reset(double _better_than = 1e20);

      // Finds the smallest strain tensor (in terms of Frobenius norm) that deforms (*this) into a lattice symmetrically equivalent to 'child_lattice'
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

      const DMatType &parent_matrix() const {
        return m_parent;
      }

      const DMatType &child_matrix() const {
        return m_child;
      }

    private:
      DMatType m_parent, m_child;
      //Conversion matrices:
      //  m_N = m_U * m_inv_count().inverse() * m_V_inv
      DMatType m_U, m_V_inv;

      StrainCostCalculator m_calc;

      // m_scale = (det(m_child)/det(m_parent))^(2/3) = det(m_F)^(2/3)
      double m_vol_factor;

      // m_atomic_vol = (det(m_child)/num_atoms)
      double m_atomic_vol;
      double m_tol;
      int m_range;

      // pointer to static list of unimodular matrices
      std::vector<Eigen::Matrix3i> const *m_mvec_ptr;

      // point group matrices, in fractional coordinates
      std::vector<Eigen::Matrix3i> m_fsym_mats;

      mutable double m_cost;
      mutable Index m_currmat;
      mutable DMatType m_F, m_N, m_dcache;
      mutable IMatType m_icache;

      ///\brief Returns the inverse of the current transformation matrix under consideration
      // We treat the unimodular matrices as the inverse of the transformation
      // matrices that we are actively considering, allowing fewer matrix inversions
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
}
#endif
