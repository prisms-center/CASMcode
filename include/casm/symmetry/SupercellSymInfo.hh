#ifndef CASM_SupercellSymInfo
#define CASM_SupercellSymInfo

#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SymGroupRepID.hh"
#include "casm/basis_set/DoFDecl.hh"
#include "casm/crystallography/PrimGrid.hh"

namespace CASM {
  class Structure;
  class PermuteIterator;

  class SupercellSymInfo {
  public:
    using permute_const_iterator = PermuteIterator;
    using SublatSymReps = std::vector<SymGroupRep::RemoteHandle>;

    SupercellSymInfo(Lattice const &_prim_lat,
                     Lattice const &_super_lat,
                     Index NB,
                     SymGroup const &_prim_factor_group,
                     SymGroupRepID basis_permutation_symrep_ID,
                     std::map<DoFKey, SymGroupRepID> const &global_dof_symrep_IDs,
                     std::vector<SymGroupRepID> const &occ_symrep_IDs,
                     std::map<DoFKey, std::vector<SymGroupRepID> > const &local_dof_symrep_IDs);

    SymGroupRep::RemoteHandle const &basis_permutation_symrep() const {
      return m_basis_perm_symrep;
    }

    // the permutation_symrep is the SymGroupRep of prim().factor_group() that describes how
    // operations of m_factor_group permute sites of the Supercell.
    SymGroupRep::RemoteHandle const &site_permutation_symrep() const;

    SymGroupRep::RemoteHandle const &global_dof_symrep(DoFKey const &_key) const {
      return m_global_dof_symreps.at(_key);
    }

    SublatSymReps const &occ_symreps() const {
      return m_occ_symreps;
    }

    SublatSymReps const &local_dof_symreps(DoFKey const &_key) const {
      return m_local_dof_symreps.at(_key);
    }

    PrimGrid const &prim_grid() const {
      return m_prim_grid;
    }

    SymGroup const &factor_group() const {
      return m_factor_group;
    }

    /// \brief Begin iterator over pure translational permutations
    permute_const_iterator translate_begin() const;

    /// \brief End iterator over pure translational permutations
    permute_const_iterator translate_end() const;

    const Permutation &factor_group_permute(Index i) const;

    // begin and end iterators for iterating over translation and factor group permutations
    permute_const_iterator permute_begin() const;
    permute_const_iterator permute_end() const;
    permute_const_iterator permute_it(Index fg_index, Index trans_index) const;
    permute_const_iterator permute_it(Index fg_index, UnitCell trans) const;

  private:
    PrimGrid m_prim_grid;

    // m_factor_group is factor group of the super cell, found by identifying the subgroup of
    // (*this).prim().factor_group() that leaves the supercell lattice vectors unchanged
    // if (*this).prim() is actually primitive, then m_factor_group.size() <= 48
    // NOTE: This is different from the SymGroup found by doing (*this).superstruc().factor_group()
    //       if Tprim is the translation group formed by the primitive cell lattice vectors, then
    //       m_factor_group is the group formed by the cosets of Tprim in the supercell space group
    //       if Tsuper is the translation group formed by the supercell lattice vectors, then,
    //       m_occupation(init_config.occupation()),
    //       m_displacement(init_config.displacement()),
    //       m_strain(init_config.supercell().strain
    //       (*this).superstruc().factor_group() is the group formed by the cosets of Tsuper in the supercell space group
    mutable SymGroup m_factor_group;

    SymGroupRep::RemoteHandle m_basis_perm_symrep;

    SublatSymReps m_occ_symreps;

    std::map<DoFKey, SublatSymReps> m_local_dof_symreps;

    std::map<DoFKey, SymGroupRep::RemoteHandle> m_global_dof_symreps;

    // m_perm_symrep_ID is the ID of the SymGroupRep of prim().factor_group() that describes how
    // operations of m_factor_group permute sites of the Supercell.
    // NOTE: The permutation representation is for (*this).prim().factor_group(), which may contain
    //       more operations than m_factor_group, so the Permutation SymGroupRep may have 'gaps' at the
    //       operations that aren't in m_factor_group. You should access elements of the SymGroupRep using
    //       the the Supercel::factor_group_permute(int) method, so that you don't encounter the gaps
    //       OR, see note for Supercell::permutation_symrep() below.
    mutable SymGroupRep::RemoteHandle m_site_perm_symrep;

  };

  template<typename IterType>
  std::vector<PermuteIterator> scel_subset_group(IterType begin, IterType end, SupercellSymInfo const &_syminfo);

  template<typename IterType>
  std::pair<MasterSymGroup, SymGroupRepID> collective_dof_symrep(IterType begin,
                                                                 IterType end,
                                                                 SupercellSymInfo const &_syminfo,
                                                                 DoFKey const &_key,
                                                                 std::vector<PermuteIterator> const &_group);

  template<typename IterType>
  Eigen::MatrixXd collective_dof_normal_coords(IterType begin,
                                               IterType end,
                                               SupercellSymInfo const &_syminfo,
                                               DoFKey const &_key,
                                               std::vector<PermuteIterator> const &_group);

}

#endif
