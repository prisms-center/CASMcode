#include "casm/basis_set/DoFIsEquivalent.hh"

#include "casm/basis_set/DoF.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/crystallography/SymType.hh"

namespace CASM {

bool DoFIsEquivalent::operator()(DoFSet const &other) const {
  return _label_equiv(other) && _vector_equiv(other.basis(), m_dof.basis());
}

bool DoFIsEquivalent::operator()(xtal::SymOp const &_op) const {
  return _vector_equiv(
      m_dof.basis(),
      DoF::BasicTraits(m_dof.type_name())
              .symop_to_matrix(get_matrix(_op), get_translation(_op),
                               get_time_reversal(_op)) *
          m_dof.basis());
}

bool DoFIsEquivalent::operator()(xtal::SymOp const &_op,
                                 DoFSet const &other) const {
  return _label_equiv(other) &&
         _vector_equiv(
             other.basis(),
             DoF::BasicTraits(m_dof.type_name())
                     .symop_to_matrix(get_matrix(_op), get_translation(_op),
                                      get_time_reversal(_op)) *
                 m_dof.basis());
}

bool DoFIsEquivalent::_label_equiv(DoFSet const &other) const {
  return m_dof.type_name() == m_dof.type_name() && m_dof.size() == other.size();
}

bool DoFIsEquivalent::_vector_equiv(
    Eigen::Ref<const Eigen::MatrixXd> const &_before_basis,
    Eigen::Ref<const Eigen::MatrixXd> const &_after_basis) const {
  if (_before_basis.cols() != _after_basis.cols()) return false;

  // If not a square matrix, make sure that column-space _after_basis is
  // identical to that of _before_basis
  if (_before_basis.cols() < _before_basis.rows()) {
    // Find rank of augmented matrix. If it is the same as rank of DoF Basis,
    // then the two matrices are similar (and, thus, equivalent)
    Eigen::MatrixXd aug(_before_basis.rows(), 2 * _before_basis.cols());
    aug << _before_basis, _after_basis;
    // std::cout << "Basis before: \n" << _before_basis << "\n\n"
    //        << "Basis after: \n" << _after_basis << "\n\n"
    //        << "DoF size: " << m_dof.size() << "\n\n"
    //        << "Augmented rank: " << aug.colPivHouseholderQr().rank() <<
    //        "\n\n";
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(aug);
    qr.setThreshold(m_tol);
    if (qr.rank() != m_dof.size()) return false;
  }

  // Find the transform m_U where
  // _before_basis*m_U = _after_basis
  m_U = _before_basis.colPivHouseholderQr().solve(_after_basis);
  if (!(m_U.transpose() * m_U).eval().isIdentity(m_tol)) {
    throw std::runtime_error(
        "Cannot find orthogonal symmetry representation for DoF \"" +
        m_dof.type_name() + "\". Please review inputs.\n");
  }
  return true;
}

}  // namespace CASM
