

namespace CASM {


  /// returns true if m_dof and _other have matching labels, and m_dof = P.permute(_other)
  template<typename OccType>
  bool OccupantDoFIsEquivalent<OccType>::operator()(OccDoFType const &_other) const {
    if(_other.size() != m_dof.size())
      return false;
    Index j;
    for(Index i = 0; i < m_dof.size(); ++i) {
      for(j = 0; j < _other.size(); ++j) {
        if(identical(m_dof[i], _other[j], m_tol)) {
          m_P[i] = j;
          break;
        }
      }
      if(j == _other.size()) {
        return false;
      }
    }
    return true;
  }

  /// returns true if copy_apply(_op,m_dof) = P.permute(m_dof)
  template<typename OccType>
  bool OccupantDoFIsEquivalent<OccType>::operator()(SymOp const &_op) const {
    Index j;
    for(Index i = 0; i < m_dof.size(); ++i) {
      OccType t_occ = copy_apply(_op, m_dof[i]);
      for(j = 0; j < m_dof.size(); ++j) {
        if(identical(t_occ, m_dof[j], m_tol)) {
          m_P[i] = j;
          break;
        }
      }
      if(j == m_dof.size()) {
        return false;
      }
    }
    return true;
  }

  /// returns true if copy_apply(_op,m_dof) =  P.permute(_other)
  template<typename OccType>
  bool OccupantDoFIsEquivalent<OccType>::operator()(SymOp const &_op, OccDoFType const &_other) const {
    if(_other.size() != m_dof.size())
      return false;
    Index j;
    for(Index i = 0; i < m_dof.size(); ++i) {
      OccType t_occ = copy_apply(_op, m_dof[i]);
      for(j = 0; j < _other.size(); ++j) {
        if(identical(t_occ, _other[j], m_tol)) {
          m_P[i] = j;
          break;
        }
      }
      if(j == _other.size()) {
        return false;
      }
    }
    return true;

  }

}
