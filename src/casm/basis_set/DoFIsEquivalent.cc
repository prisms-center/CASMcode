#include "casm/basis_set/DoFIsEquivalent.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {

  bool DoFIsEquivalent::operator()(DoFSet const &other) const {
    return _label_equiv(other) && _vector_equiv(other.basis());
  }

  bool DoFIsEquivalent::operator()(SymOp const &_op) const {

    return _vector_equiv(DoFType::traits(m_dof.type_name()).symop_to_matrix(_op) * m_dof.basis());

  }

  bool DoFIsEquivalent::operator()(SymOp const &_op, DoFSet const &other) const {
    return _label_equiv(other) && _vector_equiv(DoFType::traits(m_dof.type_name()).symop_to_matrix(_op) * other.basis());

  }

  bool DoFIsEquivalent::_label_equiv(DoFSet const &other) const {
    if(m_dof.type_name() != other.type_name() || m_dof.size() != other.size())
      return false;
    for(Index i = 0; i < m_dof.size(); i++) {
      if(m_dof[i].var_name() != other[i].var_name())
        return false;
    }
    return true;
  }

  bool DoFIsEquivalent::_vector_equiv(Eigen::Ref<const Eigen::MatrixXd> const &_basis) const {
    if(m_dof.basis().cols() != _basis.cols())
      return false;
    if(m_dof.basis().rows() == m_dof.basis().cols()) {

      m_U = _basis.colPivHouseholderQr().solve(m_dof.basis());
      if(!(m_U.transpose()*m_U).eval().isIdentity(1e-5)) {
        throw std::runtime_error("Cannot find orthogonal symmetry representation for DoF \"" + m_dof.type_name() + "\". Please review inputs.\n");
      }
      return true;
    }

    //Find rank of augmented matrix. If it is the same as rank of DoF Basis, then
    //the two matrices are similar (and, thus, equivalent)
    Eigen::MatrixXd aug(m_dof.basis().rows(), 2 * m_dof.basis().cols());
    aug << m_dof.basis(), _basis;
    //std::cout << "Basis before: \n" << m_dof.basis() << "\n\n"
    //        << "Basis after: \n" << _basis << "\n\n"
    //        << "DoF size: " << m_dof.size() << "\n\n"
    //        << "Augmented rank: " << aug.colPivHouseholderQr().rank() << "\n\n";

    if(aug.colPivHouseholderQr().rank() != m_dof.size())
      return false;
    m_U = _basis.colPivHouseholderQr().solve(m_dof.basis());

    return true;
  }

}
