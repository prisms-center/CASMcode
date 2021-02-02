#ifndef LATTICEMAP_HH
#define LATTICEMAP_HH

#include "casm/container/Counter.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/global/definitions.hh"

namespace CASM {
namespace xtal {

struct SymOp;
typedef std::vector<SymOp> SymOpVector;

/** \ingroup Lattice
 *  @{
 */

class StrainCostCalculator {
 public:
  StrainCostCalculator(Eigen::Ref<const Eigen::MatrixXd> const
                           &strain_gram_mat = Eigen::MatrixXd::Identity(9, 9));

  //\brief Isotropic strain cost, without gram matrix
  static double isotropic_strain_cost(
      Eigen::Matrix3d const &_deformation_gradient);

  //\brief Isotropic strain cost, without gram matrix
  static double isotropic_strain_cost(
      Eigen::Matrix3d const &_deformation_gradient, double _vol_factor);

  // \brief Volumetric factor :
  // pow(abs(_deformation_gradient.determinant()),1./3.), used to normalize the
  // strain cost to make it volume-independent
  static double vol_factor(Eigen::Matrix3d const &_deformation_gradient) {
    return pow(std::abs(_deformation_gradient.determinant()), 1. / 3.);
  }

  //\brief Anisotropic strain cost; utilizes stored gram matrix to compute
  // strain cost
  double strain_cost(Eigen::Matrix3d const &_deformation_gradient) const;

  //\brief Anisotropic strain cost; utilizes stored gram matrix to compute
  // strain cost
  double strain_cost(Eigen::Matrix3d const &_deformation_gradient,
                     double _vol_factor) const;

  //\brief Symmetrized strain cost; Utilizes the parent point group symmetry to
  // calculate only the symmetry breaking lattice cost
  double strain_cost(Eigen::Matrix3d const &_deformation_gradient,
                     SymOpVector const &parent_point_group) const;

 private:
  Eigen::MatrixXd m_gram_mat;
  bool m_sym_cost;

  mutable Eigen::Matrix3d m_cache;
  mutable Eigen::Matrix3d m_cache_inv;
};

/// Find the parent mapping of Lattice _parent onto Lattice _child
/// Denoting _parent.lat_column_mat() as 'parent' and _child.lat_column_mat() as
/// 'child', we want a mapping
///            child = deformation_gradient*parent*N
/// where deformation_gradient is an arbitrary 3x3 matrix and 'N' is a 3x3
/// unimodular (i.e., determinant=+/-1) integer matrix such that
/// cost(deformation_gradient) is minimized with respect to the matrix 'N' For
/// the cost function we use:
///         cost(deformation_gradient) = (_child.volume()/num_atoms)^(2/3) *
///         trace(D.transpose()*D) / g
/// where 'D' is the isovolumetric strain
///         D = deformation_gradient/det(deformation_gradient)^(1/3)-Identity
/// and 'g' is a const geometric factor.  The cost function corresponds to the
/// mean-square displacement of a point on the surface of a sphere having
/// V=_child.volume()/num_atoms (i.e., the atomic volume of the child crystal)
/// when the sphere is deformed at constant volume by
/// deformation_gradient/det(deformation_gradient)^(1/3)

class LatticeMap {
 public:
  typedef Eigen::Matrix<double, 3, 3, Eigen::DontAlign> DMatType;
  typedef Eigen::Matrix<int, 3, 3, Eigen::DontAlign> IMatType;

  LatticeMap(Lattice const &_parent, Lattice const &_child, Index _num_atoms,
             int _range, SymOpVector const &_parent_point_group,
             SymOpVector const &_child_point_group,
             Eigen::Ref<const Eigen::MatrixXd> const &strain_gram_mat =
                 Eigen::MatrixXd::Identity(9, 9),
             double _init_better_than = 1e20,
             bool _symmetrize_strain_cost = false, double _xtal_tol = TOL);

  LatticeMap(Eigen::Ref<const DMatType> const &_parent,
             Eigen::Ref<const DMatType> const &_child, Index _num_atoms,
             int _range, SymOpVector const &_parent_point_group,
             SymOpVector const &_child_point_group,
             Eigen::Ref<const Eigen::MatrixXd> const &strain_gram_mat =
                 Eigen::MatrixXd::Identity(9, 9),
             double _init_better_than = 1e20,
             bool _symmetrize_strain_cost = false, double _xtal_tol = TOL);

  void reset(double _better_than = 1e20);

  // Finds the smallest strain tensor (in terms of Frobenius norm) that deforms
  // (*this) into a lattice symmetrically equivalent to 'child_lattice'
  LatticeMap const &best_strain_mapping() const;

  LatticeMap const &next_mapping_better_than(double max_cost) const;

  double strain_cost() const { return m_cost; }

  double calc_strain_cost(const Eigen::Matrix3d &deformation_gradient) const;

  const DMatType &matrixN() const { return m_N; }

  const DMatType &deformation_gradient() const {
    return m_deformation_gradient;
  }

  const DMatType &parent_matrix() const { return m_parent; }

  const DMatType &child_matrix() const { return m_child; }

  double xtal_tol() const { return m_xtal_tol; }
  bool symmetrize_strain_cost() const { return m_symmetrize_strain_cost; }

 private:
  DMatType m_parent, m_child;
  // Conversion matrices:
  //  m_N = m_U * m_inv_count().inverse() * m_V_inv
  DMatType m_U, m_V_inv;

  StrainCostCalculator m_calc;

  // m_scale = (det(m_child)/det(m_parent))^(2/3) =
  // det(m_deformation_gradient)^(2/3)
  double m_vol_factor;

  int m_range;

  // pointer to static list of unimodular matrices
  std::vector<Eigen::Matrix3i> const *m_mvec_ptr;

  // parent point group matrices, in fractional coordinates
  std::vector<Eigen::Matrix3i> m_parent_fsym_mats;

  // parent point group in cartesian coordinates
  SymOpVector m_parent_point_group;

  // child point group matrices, in fractional coordinates
  std::vector<Eigen::Matrix3i> m_child_fsym_mats;

  // flag indicating if the symmetrized strain cost should be used while
  // searching for the best lattice maps
  bool m_symmetrize_strain_cost;
  double m_xtal_tol;

  mutable double m_cost;
  mutable Index m_currmat;
  mutable DMatType m_deformation_gradient, m_N, m_dcache;
  mutable IMatType m_icache;

  ///\brief Returns the inverse of the current transformation matrix under
  /// consideration
  // We treat the unimodular matrices as the inverse of the transformation
  // matrices that we are actively considering, allowing fewer matrix inversions
  Eigen::Matrix3i const &inv_mat() const { return (*m_mvec_ptr)[m_currmat]; }

  ///\brief Number of unimodular matrices
  Index n_mat() const { return m_mvec_ptr->size(); }

  /// \brief Returns true if current transformation is the canonical equivalent
  bool _check_canonical() const;

  LatticeMap const &_next_mapping_better_than(double max_cost) const;

  // use m_deformation_gradient to calculate strain cost
  double _calc_strain_cost() const;
};

/** @} */
}  // namespace xtal
}  // namespace CASM
#endif
