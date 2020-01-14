#ifndef CASM_SupercellSymInfo
#define CASM_SupercellSymInfo

#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/global/eigen.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SymGroupRepID.hh"
#include "casm/basis_set/DoFDecl.hh"
#include <vector>

namespace CASM {
  namespace xtal {
    class Structure;
    class UnitCell;
  }
  using xtal::Structure;
  using xtal::UnitCell;

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

    const xtal::UnitCellIndexConverter &unitcell_index_converter() const {
      return m_unitcell_to_index_converter;
    }

    const xtal::UnitCellCoordIndexConverter &unitcellcoord_index_converter() const {
      return m_unitcellcoord_to_index_converter;
    }

    const std::vector<Permutation> &translation_permutations() const {
      return m_translation_permutations;
    }

    SymGroup const &factor_group() const {
      return m_factor_group;
    }

    bool has_aniso_occs() const {
      return m_has_aniso_occs;
    }

    bool has_occupation_dofs() const {
      return m_has_occupation_dofs;
    }

    const xtal::Lattice &supercell_lattice() const {
      return m_supercell_superlattice.superlattice();
    }

    const xtal::Lattice &prim_lattice() const {
      return m_supercell_superlattice.prim_lattice();
    }

    const xtal::Superlattice &superlattice() const {
      return m_supercell_superlattice;
    }

    Eigen::Matrix3l transformation_matrix() const {
      return this->superlattice().transformation_matrix();
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

    /// Couples the primitive lattice to the supercell lattice, and knows the transformation matrix
    xtal::Superlattice m_supercell_superlattice;

    //TODO: I don't think this belongs in SupercellSymInfo, but neither did PrimGrid.
    //I'm keeping the functionality where I found it for now, but we should consider moving it elsewhere
    /// Converts between ijk (UnitCell) values and their corresponding index in an unrolled vector
    xtal::UnitCellIndexConverter m_unitcell_to_index_converter;

    //TODO: See TODO comment for m_unitcell_to_index_converter
    /// Converts between bijk (UnitCellCoord) values and their corresponding linear index
    xtal::UnitCellCoordIndexConverter m_unitcellcoord_to_index_converter;

    /// Stores the permutations associated with making translations from a lattice point to the origin
    std::vector<Permutation> m_translation_permutations;

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

    // true if there species with non-identity symreps
    bool m_has_aniso_occs;

    // true if any site has occupation DoFs
    bool m_has_occupation_dofs;

    // m_perm_symrep_ID is the ID of the SymGroupRep of prim().factor_group() that describes how
    // operations of m_factor_group permute sites of the Supercell.
    // NOTE: The permutation representation is for (*this).prim().factor_group(), which may contain
    //       more operations than m_factor_group, so the Permutation SymGroupRep may have 'gaps' at the
    //       operations that aren't in m_factor_group. You should access elements of the SymGroupRep using
    //       the the Supercel::factor_group_permute(int) method, so that you don't encounter the gaps
    //       OR, see note for Supercell::permutation_symrep() below.
    mutable SymGroupRep::RemoteHandle m_site_perm_symrep;

  };

  /// \brief Find the subgroup of supercell factor group (specified by _syminfo) that leaves a subset of sites invariant
  /// An operation is included if no site indices in the subset are mapped onto indices outside the subset after its application
  /// @param begin,end Iterator pair to list of site indices that define invariant subset
  /// @param _syminfo SupercellSymInfo object that defines all symmetry properties of supercell
  /// \result vector of PermuteIterator, referenced to _syminfo, that specify the factor group operations
  template<typename IterType>
  std::vector<PermuteIterator> scel_subset_group(IterType begin, IterType end, SupercellSymInfo const &_syminfo);


  /// \brief Find the symmetry representation for group '_group' describing the transformation of DoF '_key' among a subset of sites
  /// @param begin,end Iterator pair to list of site indices that define subset of sites of interest
  /// @param _syminfo SupercellSymInfo object that defines all symmetry properties of supercell
  /// @param _key DoFKey specifying which local DoF is of interest
  /// @param _group vector of PermuteIterators forming the group that is to be represented (this may be larger than a crystallographic factor group)
  /// \result pair containing a MasterSymGroup instantiation of _group and a SymGroupRepID that can be used to access the 'collective_dof_symrep' within
  /// the returned MasterSymGroup
  template<typename IterType>
  std::pair<MasterSymGroup, SymGroupRepID> collective_dof_symrep(IterType begin,
                                                                 IterType end,
                                                                 SupercellSymInfo const &_syminfo,
                                                                 DoFKey const &_key,
                                                                 std::vector<PermuteIterator> const &_group);

  /// \brief Find symmetry-adapted normal coordinate basis vectors for action of '_group' acting on local DoF '_key' at site indices [begin,end]
  /// @param begin,end Iterator pair to list of site indices that define subset of sites of interest
  /// @param _syminfo SupercellSymInfo object that defines all symmetry properties of supercell
  /// @param _key DoFKey specifying which local DoF is of interest
  /// @param _group vector of PermuteIterators forming the group that is to be represented (this may be larger than a crystallographic factor group)
  /// \result Orthogonal matrix comprising normal coordinate basis vectors as its rows.
  template<typename IterType>
  Eigen::MatrixXd collective_dof_normal_coords(IterType begin,
                                               IterType end,
                                               SupercellSymInfo const &_syminfo,
                                               DoFKey const &_key,
                                               std::vector<PermuteIterator> const &_group);

  /// \brief Find symmetry-adapted normal coordinate basis vectors for action of '_group' acting on local DoF '_key' at site indices [begin,end]
  /// In addition to returning the normal coordinate transformation matrix, it also returns a list of dimensions of the corresponding
  /// irreducible invariant subspaces in which the normal coordinate basis vectors reside
  /// @param begin,end Iterator pair to list of site indices that define subset of sites of interest
  /// @param _syminfo SupercellSymInfo object that defines all symmetry properties of supercell
  /// @param _key DoFKey specifying which local DoF is of interest
  /// @param _group vector of PermuteIterators forming the group that is to be represented (this may be larger than a crystallographic factor group)
  /// \result Pair containin an orthogonal matrix comprising normal coordinate basis vectors as its rows and a list of subspace dimensions, corresponding to
  /// irreducible invariant subspaces in which the normal coordinate basis vectors reside
  template<typename IterType>
  std::pair<Eigen::MatrixXd, std::vector<Index>> collective_dof_normal_coords_and_irrep_dims(IterType begin,
                                              IterType end,
                                              SupercellSymInfo const &_syminfo,
                                              DoFKey const &_key,
                                              std::vector<PermuteIterator> const &_group);

  /// \brief Find symmetry-adapted normal coordinate basis vectors for action of '_group' acting on local DoF '_key' at site indices [begin,end]
  /// In addition to returning the normal coordinate transformation matrix, it also returns a list of dimensions of the corresponding
  /// irreducible invariant subspaces in which the normal coordinate basis vectors reside
  /// @param begin,end Iterator pair to list of site indices that define subset of sites of interest
  /// @param _syminfo SupercellSymInfo object that defines all symmetry properties of supercell
  /// @param _key DoFKey specifying which local DoF is of interest
  /// @param _group vector of PermuteIterators forming the group that is to be represented (this may be larger than a crystallographic factor group)
  /// @param __subspace matrix whose columns span a subspace of allowed DoF values for selected sites
  /// \result Pair containin an orthogonal matrix comprising normal coordinate basis vectors as its rows and a list of subspace dimensions, corresponding to
  /// irreducible invariant subspaces in which the normal coordinate basis vectors reside
  template<typename IterType>
  std::pair<Eigen::MatrixXd, std::vector<Index>> collective_dof_normal_coords_and_irrep_dims(IterType begin,
                                              IterType end,
                                              SupercellSymInfo const &_syminfo,
                                              DoFKey const &_key,
                                              std::vector<PermuteIterator> const &_group,
                                              Eigen::Ref<const Eigen::MatrixXd> const &_subspace);


}

#endif
