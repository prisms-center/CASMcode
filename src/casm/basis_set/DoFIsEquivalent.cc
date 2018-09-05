#include "casm/basis_set/DoF.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM{
  bool DoFIsEquivalent::operator()(DoFSet const& other) const{
    return _label_equiv(other)&&_vector_equiv(other.basis());
  }
  
  bool DoFIsEquivalent::operator()(SymOp const& _op) const{
    return _vector_equiv(DoFType::traits(m_dof().type_name()).symop_to_matrix(_op)*m_dof.basis());

  }
  
  bool DoFIsEquivalent::operator()(SymOp const& _op, DoFSet const& other) const{
    return _label_equiv(other) && _vector_equiv(DoFType::traits(m_dof().type_name()).symop_to_matrix(_op)*other.basis());

  }

  bool DoFIsEquivalent::_label_equiv()(DoFSet const& other) const{
    if(m_dof.type_name()!=other.type_name() || m_dof.size()=other.size())
      return false;
    for(Index i=0; i<m_dof.size(); i++){
      if(m_dof[i].var_name()!=other.var_name())
        return false;
    }
    return true;
  }
  
  bool DoFIsEquivalent::_vector_equiv()(Eigen::Ref<Eigen::MatrixXd> const& _basis) const{
    if((_basis.transpose()*m_dof.basis()).colPivHouseholderQr().rank()!=m_dof.size())
      return false;
    
    m_U=qr.colPivHouseholderQr().solve(m_dof.basis())
      return true;
  }


  
}
