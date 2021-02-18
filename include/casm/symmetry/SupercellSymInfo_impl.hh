#include "casm/misc/CASM_math.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymMatrixXd.hh"
#include "casm/symmetry/SymRepTools.hh"

#ifndef CASM_SupercellSymInfo_impl
#define CASM_SupercellSymInfo_impl

namespace CASM {
template <typename IterType>
std::vector<PermuteIterator> scel_subset_group(
    IterType begin, IterType end, SupercellSymInfo const &_syminfo) {
  std::vector<PermuteIterator> result;
  std::set<Index> subset(begin, end);
  PermuteIterator perm_it = _syminfo.permute_begin(),
                  perm_end = _syminfo.permute_end();

  for (; perm_it != perm_end; ++perm_it) {
    bool add_it = true;

    for (IterType it = begin; it != end; ++it) {
      if (subset.count(perm_it.permute_ind(*it)) == 0) {
        add_it = false;
        break;
      }
    }
    if (add_it) result.push_back(perm_it);
  }
  return result;
}

//***********************************************************

template <typename IterType>
std::pair<MasterSymGroup, SymGroupRepID> collective_dof_symrep(
    IterType begin, IterType end, SupercellSymInfo const &_syminfo,
    DoFKey const &_key, std::vector<PermuteIterator> const &_group) {
  std::pair<MasterSymGroup, SymGroupRepID> result;
  if (_group.empty())
    throw std::runtime_error("Empty group passed to collective_dof_symrep()");
  result.first.set_lattice(_syminfo.supercell_lattice());
  for (PermuteIterator const &perm : _group) {
    result.first.push_back(perm.sym_op());
  }

  result.second = result.first.allocate_representation();
  SupercellSymInfo::SublatSymReps const &subreps =
      _key == "occ" ? _syminfo.occ_symreps() : _syminfo.local_dof_symreps(_key);

  // make map of site_index -> beginning row in basis for that site
  // (number of rows per site == dof dimension on that site)
  std::map<Index, Index> site_index_to_basis_index;
  Index total_dim = 0;
  for (IterType it = begin; it != end; ++it) {
    Index b = _syminfo.unitcellcoord_index_converter()(*it).sublattice();
    Index site_dof_dim = subreps[b].dim();
    site_index_to_basis_index[*it] = total_dim;
    total_dim += site_dof_dim;
  }

  // make matrix rep, by filling in blocks with site dof symreps
  Eigen::MatrixXd trep(total_dim, total_dim);
  Index g = 0;
  for (PermuteIterator const &perm : _group) {
    trep.setZero();
    for (IterType it = begin; it != end; ++it) {
      Index b = _syminfo.unitcellcoord_index_converter()(*it).sublattice();
      auto ptr = (subreps[b][perm.factor_group_index()]->MatrixXd());

      // can't fail, because it was built from [begin, end)
      Index row = site_index_to_basis_index.find(*it)->second;

      // could fail, if mismatch between [begin, end) and group
      auto col_it = site_index_to_basis_index.find(perm.permute_ind(*it));
      if (col_it == site_index_to_basis_index.end()) {
        throw std::runtime_error(
            "Error in collective_dof_symrep: Input group includes permutations "
            "between selected and unselected sites.");
      }
      Index col = col_it->second;

      trep.block(row, col, ptr->rows(), ptr->cols()) = *ptr;
    }
    result.first[g++].set_rep(result.second, SymMatrixXd(trep));
  }
  result.first.sort();
  return result;
}
//***********************************************************

template <typename IterType>
std::vector<SymRepTools::IrrepInfo> irrep_decomposition(
    IterType begin, IterType end, SupercellSymInfo const &_syminfo,
    DoFKey const &_key, std::vector<PermuteIterator> const &_group) {
  Index dim =
      _syminfo.local_dof_symreps(_key)[0].dim() * std::distance(begin, end);
  return irrep_decomposition(begin, end, _syminfo, _key, _group,
                             Eigen::MatrixXd::Identity(dim, dim));
}

template <typename IterType>
std::vector<SymRepTools::IrrepInfo> irrep_decomposition(
    IterType begin, IterType end, SupercellSymInfo const &_syminfo,
    DoFKey const &_key, std::vector<PermuteIterator> const &_group,
    Eigen::Ref<const Eigen::MatrixXd> const &_subspace) {
  std::pair<MasterSymGroup, SymGroupRepID> rep_info =
      collective_dof_symrep(begin, end, _syminfo, _key, _group);

  SymGroupRep const &rep = rep_info.first.representation(rep_info.second);

  return irrep_decomposition(rep, rep_info.first, _subspace, false);
}

}  // namespace CASM

#endif
