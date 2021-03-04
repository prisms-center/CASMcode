#include "casm/misc/CASM_math.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymMatrixXd.hh"
#include "casm/symmetry/SymRepTools.hh"
#include "casm/symmetry/io/stream/SymInfo_stream_io.hh"

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
  // To build the collective DoF symrep matrices, we need to know the
  // conventions for permutations among sites, and the conventions for storing
  // site DoF symmetry representations.
  //
  // For permutation among sites, by convention:
  //     after[i] = before[perm.permute_ind(i)],
  // where:
  // - perm is a PermuteIterator
  //
  // For transforming site DoF values, by convention:
  //     op.matrix() * B_from = B_to * U(from_b, op.index()),
  // where:
  // - op.matrix(), factor group operation symmetry matrix
  // - B_from: site dof basis on "before" site
  // - B_to: site dof basis on "after" site
  // - U(from_b, op.index()): symmetry representation matrix, stored in
  //   `_syminfo.local_dof_symreps(_key)[from_b][op.index()]`
  //
  // Relationships between the site DoF on the "from" site before symmetry
  // application, to the site DoF value on the "to" site after symmetry
  // application:
  //
  //     v_standard_after = op.matrix() * v_standard_before
  //                      = op.matrix() * B_from * v_prim_from_before
  //                      = B_to * v_prim_to_after
  //     v_prim_to_after = B_to^-1 * op.matrix() * B_from * v_prim_from_before
  //     v_prim_to_after = U(from_b, op.index()) * v_prim_from_before

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
    std::cout << "---" << std::endl;
    std::cout << "building " << g << std::endl;
    std::cout << "symop: "
              << brief_description(perm.sym_op(), _syminfo.prim_lattice(),
                                   SymInfoOptions{CART})
              << std::endl;
    trep.setZero();
    for (IterType it = begin; it != end; ++it) {
      // can't fail, because it was built from [begin, end)
      Index to_site_index = *it;
      Index row = site_index_to_basis_index.find(to_site_index)->second;
      std::cout << "row: " << row << std::endl;

      // could fail, if mismatch between [begin, end) and group
      Index from_site_index = perm.permute_ind(*it);
      auto col_it = site_index_to_basis_index.find(from_site_index);
      if (col_it == site_index_to_basis_index.end()) {
        throw std::runtime_error(
            "Error in collective_dof_symrep: Input group includes permutations "
            "between selected and unselected sites.");
      }
      Index col = col_it->second;
      std::cout << "col: " << col << std::endl;

      Index from_site_b =
          _syminfo.unitcellcoord_index_converter()(from_site_index)
              .sublattice();
      Eigen::MatrixXd U =
          *(subreps[from_site_b][perm.factor_group_index()]->MatrixXd());
      std::cout << "from_site_b: " << from_site_b << std::endl;
      std::cout << "perm.factor_group_index(): " << perm.factor_group_index()
                << std::endl;
      std::cout << "*it: " << *it << std::endl;
      std::cout << "perm.permute_ind(*it): " << perm.permute_ind(*it)
                << std::endl;

      std::cout << "subreps matrix: \n" << U << std::endl;
      std::cout << std::endl;
      trep.block(row, col, U.rows(), U.cols()) = U;
    }
    std::cout << "trep: \n" << trep << std::endl;
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
