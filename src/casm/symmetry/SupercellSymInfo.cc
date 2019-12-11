#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/symmetry/SupercellSymInfo.hh"
#include "casm/symmetry/SupercellSymInfo_impl.hh"
#include "casm/symmetry/SymBasisPermute.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymPermutation.hh"
#include "casm/symmetry/SymTools.hh"

namespace CASM {

  SupercellSymInfo::SupercellSymInfo(Lattice const &_prim_lat,
                                     Lattice const &_super_lat,
                                     Index num_sites_in_prim,
                                     SymGroup const &_prim_factor_group,
                                     SymGroupRepID basis_permutation_symrep_ID,
                                     std::map<DoFKey, SymGroupRepID> const &global_dof_symrep_IDs,
                                     std::vector<SymGroupRepID> const &occ_symrep_IDs,
                                     std::map<DoFKey, std::vector<SymGroupRepID> > const &local_dof_symrep_IDs) :
    m_supercell_superlattice(_prim_lat, _super_lat),
    m_unitcell_to_index_converter(m_supercell_superlattice.transformation_matrix()),
    m_unitcellcoord_to_index_converter(m_supercell_superlattice.transformation_matrix(), num_sites_in_prim),
    m_prim_grid(_prim_lat, _super_lat, num_sites_in_prim),
    m_factor_group(sym::invariant_subgroup(_prim_factor_group, _super_lat)),
    m_basis_perm_symrep(factor_group(), basis_permutation_symrep_ID),
    m_has_aniso_occs(false),
    m_has_occupation_dofs(false) {
    for(auto const &dofID : global_dof_symrep_IDs)
      m_global_dof_symreps.emplace(std::make_pair(dofID.first, SymGroupRep::RemoteHandle(factor_group(), dofID.second)));

    for(auto const &dofID : local_dof_symrep_IDs) {
      SublatSymReps treps(num_sites_in_prim);
      for(Index b = 0; b < num_sites_in_prim; ++b) {
        if(!dofID.second[b].empty())
          treps[b] = SymGroupRep::RemoteHandle(factor_group(), dofID.second[b]);
      }
      m_local_dof_symreps.emplace(std::make_pair(dofID.first, std::move(treps)));
    }

    m_occ_symreps.resize(num_sites_in_prim);
    for(Index b = 0; b < num_sites_in_prim; ++b) {
      if(!occ_symrep_IDs[b].is_identity()) {
        m_has_aniso_occs = true;
        m_has_occupation_dofs = true;
      }
      else if(occ_symrep_IDs[b].rep_index() > 1) {
        m_has_occupation_dofs = true;
      }

      m_occ_symreps[b] = SymGroupRep::RemoteHandle(factor_group(), occ_symrep_IDs[b]);
    }
  }


  SymGroupRepID make_permutation_representation(const SymGroup &group, const xtal::LinearIndexConverter bijk_index_converter, const Structure &prim) {
    SymGroupRepID perm_rep_ID = group.allocate_representation();
    long total_sites = bijk_index_converter.total_sites();
    for(Index operation_ix = 0; operation_ix < group.size(); ++operation_ix) {
      const auto &operation = group[operation_ix];
      std::vector<Index> permutation(total_sites);
      for(Index old_l = 0; old_l < total_sites; ++old_l) {
        UnitCellCoord old_ucc = bijk_index_converter[old_l];
        UnitCellCoord new_ucc = sym::copy_apply(operation, old_ucc, prim);
        Index new_l = bijk_index_converter[new_ucc];
        permutation[new_l] = old_l;
      }

      group[operation_ix].set_rep(perm_rep_ID, SymPermutation(permutation));
    }
    return perm_rep_ID;
  }


  // the permutation_symrep is the SymGroupRep of prim().factor_group() that describes how
  // operations of m_factor_group permute sites of the Supercell.
  // NOTE: The permutation representation is for (*this).prim().factor_group(), which may contain
  //       more operations than m_factor_group, so the Permutation SymGroupRep may have 'gaps' at the
  //       operations that aren't in m_factor_group. You should access elements of the SymGroupRep using
  //       SymGroupRep::get_representation(m_factor_group[i]) or SymGroupRep::get_permutation(m_factor_group[i]),
  //       so that you don't encounter the gaps (i.e., the representation can be indexed using the
  //       SymOps of m_factor_group
  SymGroupRep::RemoteHandle const &SupercellSymInfo::site_permutation_symrep() const {
    if(m_site_perm_symrep.empty()) {
      m_site_perm_symrep = SymGroupRep::RemoteHandle(factor_group(), prim_grid().make_permutation_representation(factor_group(), basis_permutation_symrep().symrep_ID()));
    }

    return m_site_perm_symrep;
  }


  // permutation_symrep() populates permutation symrep if needed
  const Permutation &SupercellSymInfo::factor_group_permute(Index i) const {
    return *(site_permutation_symrep()[i]->permutation());
  }

  // PrimGrid populates translation permutations if needed
  //const Permutation &SupercellSymInfo::translation_permute(Index i) const {
  //return prim_grid().translation_permutation(i);
  //}

  // PrimGrid populates translation permutations if needed
  //const std::vector<Permutation> &SupercellSymInfo::translation_permute() const {
  //return prim_grid().translation_permutations();
  // }

  /// \brief Begin iterator over translation permutations
  SupercellSymInfo::permute_const_iterator SupercellSymInfo::translate_begin() const {
    return permute_begin();
  }

  /// \brief End iterator over translation permutations
  SupercellSymInfo::permute_const_iterator SupercellSymInfo::translate_end() const {
    return permute_begin().begin_next_fg_op();
  }

  SupercellSymInfo::permute_const_iterator SupercellSymInfo::permute_begin() const {
    return permute_it(0, 0); // starting indices
  }

  SupercellSymInfo::permute_const_iterator SupercellSymInfo::permute_end() const {
    return permute_it(factor_group().size(), 0); // one past final indices
  }

  SupercellSymInfo::permute_const_iterator SupercellSymInfo::permute_it(Index fg_index, Index trans_index) const {
    return permute_const_iterator(*this,
                                  fg_index, trans_index);
  }

  SupercellSymInfo::permute_const_iterator SupercellSymInfo::permute_it(Index fg_index, UnitCell trans) const {
    Index trans_index = this->unitcell_index_converter()[trans];
    return permute_it(fg_index, trans_index);
  }


}
