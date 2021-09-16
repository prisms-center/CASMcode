#include "casm/crystallography/OccupantDoFIsEquivalent.hh"

#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/SymType.hh"

namespace CASM {
namespace xtal {

/// returns true if m_dof and _other have matching labels, and m_dof =
/// P.permute(_other)
bool OccupantDoFIsEquivalent::operator()(
    std::vector<xtal::Molecule> const &_other) const {
  if (_other.size() != m_dof.size()) return false;
  Index j;
  for (Index i = 0; i < m_dof.size(); ++i) {
    for (j = 0; j < _other.size(); ++j) {
      if (m_dof[i].identical(_other[j], m_tol)) {
        m_P.set(i) = j;
        break;
      }
    }
    if (j == _other.size()) {
      return false;
    }
  }
  return true;
}

/// returns true if copy_apply(_op,m_dof) = P.permute(m_dof)
bool OccupantDoFIsEquivalent::operator()(xtal::SymOp const &_op) const {
  Index j;
  for (Index i = 0; i < m_dof.size(); ++i) {
    xtal::Molecule t_occ = sym::copy_apply(_op, m_dof[i]);
    for (j = 0; j < m_dof.size(); ++j) {
      if (t_occ.identical(m_dof[j], m_tol)) {
        m_P.set(i) = j;
        break;
      }
    }
    if (j == m_dof.size()) {
      return false;
    }
  }
  return true;
}

/// returns true if copy_apply(_op,m_dof) =  P.permute(_other)
bool OccupantDoFIsEquivalent::operator()(
    xtal::SymOp const &_op, std::vector<xtal::Molecule> const &_other) const {
  if (_other.size() != m_dof.size()) return false;
  Index j;
  for (Index i = 0; i < m_dof.size(); ++i) {
    xtal::Molecule t_occ = sym::copy_apply(_op, m_dof[i]);
    for (j = 0; j < _other.size(); ++j) {
      if (t_occ.identical(_other[j], m_tol)) {
        m_P.set(i) = j;
        break;
      }
    }
    if (j == _other.size()) {
      return false;
    }
  }
  return true;
}

}  // namespace xtal
}  // namespace CASM
