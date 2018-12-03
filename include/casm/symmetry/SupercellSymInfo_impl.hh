#include "casm/symmetry/PermuteIterator.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymMatrixXd.hh"

#ifndef CASM_SupercellSymInfo_impl
#define CASM_SupercellSymInfo_impl


namespace CASM {
  template<typename IterType>
  std::vector<PermuteIterator> scel_subset_group(IterType begin, IterType end, SupercellSymInfo const &_syminfo) {
    std::vector<PermuteIterator> result;
    std::set<Index> subset(begin, end);
    PermuteIterator perm_it = _syminfo.permute_begin(), perm_end = _syminfo.permute_end();

    for(; perm_it != perm_end; ++perm_it) {
      bool add_it = true;

      for(IterType it = begin; it != end; ++it) {
        if(subset.count(perm_it.permute_ind(*it)) == 0) {
          add_it = false;
          break;
        }
      }
      if(add_it)
        result.push_back(perm_it);
    }

  }

  //***********************************************************

  template<typename IterType>
  std::pair<MasterSymGroup, SymGroupRepID> collective_dof_symrep(IterType begin,
                                                                 IterType end,
                                                                 SupercellSymInfo const &_syminfo,
                                                                 DoFKey const &_key,
                                                                 std::vector<PermuteIterator> const &_group) {


    std::pair<MasterSymGroup, SymGroupRepID> result;// = std::make_pair(make_sym_group(_group),SymGroupRepID());
    if(_group.empty())
      throw std::runtime_error("Empty group passed to collective_dof_symrep()");
    result.first.set_lattice(_group[0].prim_grid().scel_lattice());
    for(PermuteIterator const &perm : _group) {
      result.first.push_back(perm.sym_op());
    }

    result.second = result.first.allocate_representation();
    Index subdim = 0;
    SupercellSymInfo::SublatSymReps const &subreps = _key == "occ" ? _syminfo.occ_symreps() : _syminfo.local_dof_symreps(_key);

    for(auto const &rep : subreps)
      subdim = max(subdim, rep.dim());

    //std::cout << "subdim is " << subdim << "\n";

    Index Nsite = std::distance(begin, end);

    //std::cout << "Nsite is " << Nsite << "\n";

    Eigen::MatrixXd trep(subdim * Nsite, subdim * Nsite);
    Index g = 0;
    for(PermuteIterator const &perm : _group) {
      trep.setZero();
      for(IterType it = begin; it != end; ++it) {
        Index b = _syminfo.prim_grid().sublat(*it);
        auto ptr = (subreps[b][perm.factor_group_index()]->MatrixXd());
        trep.block(subdim * (*it), subdim * perm.permute_ind(*it), ptr->rows(), ptr->cols()) = *ptr;
      }
      std::cout << "trep " << g << " is \n" << trep << "\n\n";
      result.first[g++].set_rep(result.second, SymMatrixXd(trep));
      //std::cout << "cartesian symop is: \n" << perm.sym_op().matrix() << "\n";
    }
    result.first.sort();
    return result;
  }
  //***********************************************************

  template<typename IterType>
  Eigen::MatrixXd collective_dof_normal_coords(IterType begin,
                                               IterType end,
                                               SupercellSymInfo const &_syminfo,
                                               DoFKey const &_key,
                                               std::vector<PermuteIterator> const &_group) {
    return collective_dof_normal_coords_and_irrep_dims(begin, end, _syminfo, _key, _group).first;

  }

  template<typename IterType>
  std::pair<Eigen::MatrixXd, std::vector<Index>> collective_dof_normal_coords_and_irrep_dims(IterType begin,
                                              IterType end,
                                              SupercellSymInfo const &_syminfo,
                                              DoFKey const &_key,
  std::vector<PermuteIterator> const &_group) {
    //std::vector<PermuteIterator> perms = scel_subset_group(begin,end,_syminfo);
    std::pair<MasterSymGroup, SymGroupRepID> rep_info = collective_dof_symrep(begin, end, _syminfo, _key, _group);

    SymGroupRep::RemoteHandle rep(rep_info.first, rep_info.second);
    return rep->_get_irrep_trans_mat_blind(rep_info.first);

  }
}

#endif
