#include <casm/crystallography/DoFSet.hh>
#include <casm/crystallography/SymType.hh>
#include <unordered_set>
#include <vector>

#include "casm/crystallography/Site.hh"

namespace CASM {
namespace xtal {
/// Returns descriptive names of the components in a DoFSet, using
/// AnisoValTraits::variable_descriptors()
std::vector<std::string> component_descriptions(DoFSet const &dofset) {
  if (dofset.basis().isIdentity(TOL) &&
      dofset.basis().rows() == dofset.basis().cols())
    return dofset.traits().variable_descriptions();
  else
    return dofset.component_names();
}

bool DoFSetIsEquivalent_f::_traits_match(const DoFSet &other_value) const {
  return m_reference_dofset.traits().name() == other_value.traits().name();
}

bool DoFSetIsEquivalent_f::_basis_spans_same_space(
    const DoFSet &other_value) const {
  const auto &reference_basis = m_reference_dofset.basis();
  const auto &other_basis = other_value.basis();

  if (reference_basis.cols() != other_basis.cols()) return false;

  // If not a square matrix, make sure that column-space _after_basis is
  // identical to that of _before_basis
  if (reference_basis.cols() < other_basis.rows()) {
    // Find rank of augmented matrix. If it is the same as rank of DoF Basis,
    // then the two matrices are similar (and, thus, equivalent)
    Eigen::MatrixXd aug(reference_basis.rows(), 2 * reference_basis.cols());
    aug << reference_basis, other_basis;

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(aug);
    qr.setThreshold(m_tol);
    if (qr.rank() != m_reference_dofset.dim()) return false;
  }

  return true;
}

bool DoFSetIsEquivalent_f::operator()(const DoFSet &other_value) const {
  return this->_traits_match(other_value) &&
         this->_basis_spans_same_space(other_value);
}

Eigen::MatrixXd dofset_transformation_matrix(const Eigen::MatrixXd &from_basis,
                                             const Eigen::MatrixXd &to_basis,
                                             double tol) {
  Eigen::MatrixXd U = from_basis.colPivHouseholderQr().solve(to_basis);
  if (!(U.transpose() * U).eval().isIdentity(tol)) {
    throw std::runtime_error("Cannot find orthogonal symmetry representation!");
  }
  return U;
}
}  // namespace xtal
}  // namespace CASM

namespace CASM {
namespace sym {

/// \brief Copy and apply SymOp to a DoFSet
xtal::DoFSet copy_apply(const xtal::SymOp &op, const xtal::DoFSet &dof) {
  Eigen::MatrixXd transformation_matrix = dof.traits().symop_to_matrix(
      get_matrix(op), get_translation(op), get_time_reversal(op));
  return xtal::DoFSet(dof.traits(), dof.component_names(),
                      transformation_matrix * dof.basis());
}

/// \brief Copy and apply SymOp to a SiteDoFSet
xtal::SiteDoFSet copy_apply(const xtal::SymOp &op,
                            const xtal::SiteDoFSet &dof) {
  Eigen::MatrixXd transformation_matrix = dof.traits().symop_to_matrix(
      get_matrix(op), get_translation(op), get_time_reversal(op));
  return xtal::SiteDoFSet(dof.traits(), dof.component_names(),
                          transformation_matrix * dof.basis(),
                          dof.excluded_occupants());
}

}  // namespace sym

}  // namespace CASM
