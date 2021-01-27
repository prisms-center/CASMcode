#include "casm/symmetry/SymOpRepresentation.hh"

#include <memory>

#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SymMatrixXd.hh"
#include "casm/symmetry/SymPermutation.hh"

namespace CASM {
//*******************************************************************************************

std::unique_ptr<SymOpRepresentation> SymOpRepresentation::inverse() const {
  SymOpRepresentation *res = inverse_impl();
  res->m_master_group = m_master_group;
  res->m_rep_ID = m_rep_ID;
  if (has_valid_master()) res->m_op_index = ind_inverse();
  return std::unique_ptr<SymOpRepresentation>(res);
}

//*******************************************************************************************
SymOpRepresentation const &SymOpRepresentation::representation(
    SymGroupRepID _rep_ID) const {
  assert(has_valid_master() && !_rep_ID.empty());
  return *(master_group().representation(_rep_ID)[index()]);
}

//*******************************************************************************************
Eigen::MatrixXd const *SymOpRepresentation::get_matrix_rep(
    SymGroupRepID _rep_ID) const {
  assert(has_valid_master() && !_rep_ID.empty());
  return (master_group().representation(_rep_ID)[index()])->MatrixXd();
}

//*******************************************************************************************

SymBasisPermute const *SymOpRepresentation::get_basis_permute_rep(
    SymGroupRepID _rep_ID) const {
  assert(has_valid_master() && !_rep_ID.empty());
  return (master_group().representation(_rep_ID)[index()])->ucc_permutation();
}
//*******************************************************************************************

Permutation const *SymOpRepresentation::get_permutation_rep(
    SymGroupRepID _rep_ID) const {
  assert(has_valid_master() && !_rep_ID.empty());
  return (master_group().representation(_rep_ID)[index()])->permutation();
}

//*******************************************************************************************

Array<Eigen::MatrixXd const *> SymOpRepresentation::get_matrix_reps(
    Array<SymGroupRepID> _rep_IDs) const {
  Array<Eigen::MatrixXd const *> tmat;
  for (Index i = 0; i < _rep_IDs.size(); i++) {
    tmat.push_back(get_matrix_rep(_rep_IDs[i]));
  }
  return tmat;
}

//**********************************************************
void SymOpRepresentation::set_rep(SymGroupRepID _rep_ID,
                                  const SymOpRepresentation &op_rep) const {
  assert(has_valid_master() && !_rep_ID.empty());
  return master_group().set_rep(_rep_ID, op_rep, index());
}

//*******************************************************************************************

void SymOpRepresentation::set_identifiers(const MasterSymGroup &new_group,
                                          SymGroupRepID new_rep_ID) {
  m_master_group = &new_group;
  m_rep_ID = new_rep_ID;
  SymGroupRep const &trep(new_group.representation(m_rep_ID));
  Index i;
  for (i = 0; i < trep.size(); i++) {
    if (this == trep[i]) {
      m_op_index = i;
      break;
    }
  }

  if (i == new_group.size()) m_op_index = -1;

  return;
}

//*******************************************************************************************

void SymOpRepresentation::set_identifiers(const MasterSymGroup &new_group,
                                          SymGroupRepID new_rep_ID,
                                          Index new_op_index) {
  m_master_group = &new_group;
  m_rep_ID = new_rep_ID;
  m_op_index = new_op_index;

  return;
}

//*******************************************************************************************

Index SymOpRepresentation::ind_inverse() const {
  assert(
      has_valid_master() &&
      "In SymOpRepresentation::ind_inverse(), head_group is uninitialized!!");
  return master_group().ind_inverse(index());
}

//*******************************************************************************************

Index SymOpRepresentation::ind_prod(const SymOpRepresentation &RHS) const {
  assert(has_valid_master() &&
         "In SymOpRepresentation::ind_prod(), head_group is uninitialized!!");
  return master_group().ind_prod(index(), RHS.index());
}

//*******************************************************************************************

}  // namespace CASM
