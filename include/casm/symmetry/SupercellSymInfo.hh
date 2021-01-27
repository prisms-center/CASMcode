#ifndef CASM_SupercellSymInfo
#define CASM_SupercellSymInfo

#include <vector>

#include "casm/crystallography/DoFDecl.hh"
#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/crystallography/Superlattice.hh"
#include "casm/global/eigen.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SymGroupRepID.hh"

namespace CASM {
namespace SymRepTools {
struct IrrepInfo;
}

namespace xtal {
class UnitCell;
}
using xtal::UnitCell;

class PermuteIterator;

/// \brief A class that collects all symmetry information for for performing
/// symmetry transformations on the site indices, site DoFs, and global DoFs of
/// a Supercell or Configuration
///
class SupercellSymInfo {
 public:
  using permute_const_iterator = PermuteIterator;
  using SublatSymReps = std::vector<SymGroupRep::RemoteHandle>;

  ///\brief Construct with primitive and super lattice, number of sublattice and
  ///all relevant representation IDs
  SupercellSymInfo(Lattice const &_prim_lat, Lattice const &_super_lat,
                   Index number_of_sublats, SymGroup const &_prim_factor_group,
                   SymGroupRepID basis_permutation_symrep_ID,
                   std::map<DoFKey, SymGroupRepID> const &global_dof_symrep_IDs,
                   std::vector<SymGroupRepID> const &occ_symrep_IDs,
                   std::map<DoFKey, std::vector<SymGroupRepID> > const
                       &local_dof_symrep_IDs);

  /// Returns a "RemoteHandle" to the sublattice permutation representation of
  /// the supercell's factor group
  SymGroupRep::RemoteHandle const &basis_permutation_symrep() const;

  /// Returns a "RemoteHandle" to the site permutation representation of the
  /// supercell's factor group
  SymGroupRep::RemoteHandle const &site_permutation_symrep() const;

  /// Returns a "RemoteHandle" to the global DoF matrix representation of the
  /// supercell's factor group
  SymGroupRep::RemoteHandle const &global_dof_symrep(DoFKey const &_key) const {
    return m_global_dof_symreps.at(_key);
  }

  /// \brief SymGroupRep handle for occupant permutation representation of
  /// supercell's factor group An occupant permutation describes how the
  /// occupants at a site change identity due to a spatial transformation
  /// EXAMPLE: The label of a x-oriented dimer and a y-oriented dimer will be
  /// exchanged due to a 90-degree rotation about z
  ///          If these are the only two species at the site, the permutation
  ///          will have the form p=[1,0], with p[index_after]=index_before If
  ///          sites are also permuted, then the occ_symrep is used as
  ///              site_at(site_index_after).specie(index_after)=site_at(site_index_before).specie(p[index_after]);
  SublatSymReps const &occ_symreps() const { return m_occ_symreps; }

  /// \brief SymGroupRep handle for site DoF matrix representation of
  /// supercell's factor group
  ///          If sites are also permuted, then the local_dof_symrep is used as
  ///              site_at(site_index_after).dof_value(_key) =
  ///              representation_matrix *
  ///              site_at(site_index_before).dof_value(_key)
  SublatSymReps const &local_dof_symreps(DoFKey const &_key) const {
    return m_local_dof_symreps.at(_key);
  }

  /// \brief UnitCellIndexConverter for this superlattice/primlattice pair
  /// Used to convert from lattice translations to site permutations
  const xtal::UnitCellIndexConverter &unitcell_index_converter() const {
    return m_unitcell_to_index_converter;
  }

  /// \brief UnitCellCoordIndexConverter for this superstructure/primstructure
  /// pair Used to convert from lattice translations to site permutations
  const xtal::UnitCellCoordIndexConverter &unitcellcoord_index_converter()
      const {
    return m_unitcellcoord_to_index_converter;
  }

  /// \brief Permutations describing reordering of sites of supercell due to a
  /// lattice translation of the primitive translation_permutation()[l] gives
  /// the site permutation due to translating to origin cell to unitcell[l] of
  /// the supercell
  const std::vector<Permutation> &translation_permutations() const {
    return m_translation_permutations;
  }

  /// \brief Subgroup of primitive-cell factor group operations that leave
  /// supercell lattice invariant
  SymGroup const &factor_group() const { return m_factor_group; }

  /// \brief true if any species are anisotropic, such that the occ_symreps are
  /// non-trivial
  bool has_aniso_occs() const { return m_has_aniso_occs; }

  /// \brief true if any sublattice has more than one allowed occupant
  bool has_occupation_dofs() const { return m_has_occupation_dofs; }

  /// \brief const reference to supercell lattice
  const xtal::Lattice &supercell_lattice() const {
    return m_supercell_superlattice.superlattice();
  }

  /// \brief const reference to primitive lattice
  const xtal::Lattice &prim_lattice() const {
    return m_supercell_superlattice.prim_lattice();
  }

  /// \brief const reference to superlattice
  const xtal::Superlattice &superlattice() const {
    return m_supercell_superlattice;
  }

  /// \brief long-int transformation from primitive lattice vectors to supercell
  /// lattice vectors
  ///   supercell_lattice().lat_column_mat() = prim_lattice().lat_column_mat() *
  ///   transformation_matrix_to_super()
  Eigen::Matrix3l transformation_matrix_to_super() const {
    return this->superlattice().transformation_matrix_to_super();
  }

  /// \brief Begin PermuteIterator over pure translational permutations
  //// Equivalent to permute_begin()
  permute_const_iterator translate_begin() const;

  /// End PermuteIterator over pure translational permutations
  permute_const_iterator translate_end() const;

  /// Site permutation corresponding to supercell factor group operation
  const Permutation &factor_group_permute(
      Index supercell_factor_group_index) const;

  // begin and end iterators for iterating over translation and factor group
  // permutations
  permute_const_iterator permute_begin() const;
  permute_const_iterator permute_end() const;
  permute_const_iterator permute_it(Index supercell_factor_group_index,
                                    Index translation_index) const;
  permute_const_iterator permute_it(Index supercell_factor_group_index,
                                    UnitCell translation) const;

 private:
  /// Couples the primitive lattice to the supercell lattice, and knows the
  /// transformation matrix
  xtal::Superlattice m_supercell_superlattice;

  // TODO: I don't think this belongs in SupercellSymInfo, but neither did
  // PrimGrid. I'm keeping the functionality where I found it for now, but we
  // should consider moving it elsewhere
  /// Converts between ijk (UnitCell) values and their corresponding index in an
  /// unrolled vector
  xtal::UnitCellIndexConverter m_unitcell_to_index_converter;

  // TODO: See TODO comment for m_unitcell_to_index_converter
  /// Converts between bijk (UnitCellCoord) values and their corresponding
  /// linear index
  xtal::UnitCellCoordIndexConverter m_unitcellcoord_to_index_converter;

  /// Stores the permutations associated with making translations from a lattice
  /// point to the origin
  std::vector<Permutation> m_translation_permutations;

  // m_factor_group is factor group of the super cell, found by identifying the
  // subgroup of
  // (*this).prim().factor_group() that leaves the supercell lattice vectors
  // unchanged if (*this).prim() is actually primitive, then
  // m_factor_group.size() <= 48 NOTE: This is different from the SymGroup found
  // by doing (*this).superstruc().factor_group()
  //       if Tprim is the translation group formed by the primitive cell
  //       lattice vectors, then m_factor_group is the group formed by the
  //       cosets of Tprim in the supercell space group if Tsuper is the
  //       translation group formed by the supercell lattice vectors, then,
  //       m_occupation(init_config.occupation()),
  //       m_displacement(init_config.displacement()),
  //       m_strain(init_config.supercell().strain
  //       (*this).superstruc().factor_group() is the group formed by the cosets
  //       of Tsuper in the supercell space group
  mutable SymGroup m_factor_group;

  SymGroupRep::RemoteHandle m_basis_perm_symrep;

  SublatSymReps m_occ_symreps;

  std::map<DoFKey, SublatSymReps> m_local_dof_symreps;

  std::map<DoFKey, SymGroupRep::RemoteHandle> m_global_dof_symreps;

  // true if there species with non-identity symreps
  bool m_has_aniso_occs;

  // true if any site has occupation DoFs
  bool m_has_occupation_dofs;

  // m_perm_symrep_ID is the ID of the SymGroupRep of prim().factor_group() that
  // describes how operations of m_factor_group permute sites of the Supercell.
  // NOTE: The permutation representation is for (*this).prim().factor_group(),
  // which may contain
  //       more operations than m_factor_group, so the Permutation SymGroupRep
  //       may have 'gaps' at the operations that aren't in m_factor_group. You
  //       should access elements of the SymGroupRep using the the
  //       Supercel::factor_group_permute(int) method, so that you don't
  //       encounter the gaps OR, see note for Supercell::permutation_symrep()
  //       below.
  mutable SymGroupRep::RemoteHandle m_site_perm_symrep;
};

/// \brief Find the subgroup of supercell factor group (specified by _syminfo)
/// that leaves a subset of sites invariant An operation is included if no site
/// indices in the subset are mapped onto indices outside the subset after its
/// application
/// @param begin,end Iterator pair to list of site indices that define invariant
/// subset
/// @param _syminfo SupercellSymInfo object that defines all symmetry properties
/// of supercell \result vector of PermuteIterator, referenced to _syminfo, that
/// specify the factor group operations
template <typename IterType>
std::vector<PermuteIterator> scel_subset_group(
    IterType begin, IterType end, SupercellSymInfo const &_syminfo);

/// \brief Find the symmetry representation for group '_group' describing the
/// transformation of DoF '_key' among a subset of sites
/// @param begin,end Iterator pair to list of site indices that define subset of
/// sites of interest
/// @param _syminfo SupercellSymInfo object that defines all symmetry properties
/// of supercell
/// @param _key DoFKey specifying which local DoF is of interest
/// @param _group vector of PermuteIterators forming the group that is to be
/// represented (this may be larger than a crystallographic factor group)
/// \result pair containing a MasterSymGroup instantiation of _group and a
/// SymGroupRepID that can be used to access the 'collective_dof_symrep' within
/// the returned MasterSymGroup
template <typename IterType>
std::pair<MasterSymGroup, SymGroupRepID> collective_dof_symrep(
    IterType begin, IterType end, SupercellSymInfo const &_syminfo,
    DoFKey const &_key, std::vector<PermuteIterator> const &_group);

/// \brief Find symmetry-adapted normal coordinate basis vectors for action of
/// '_group' acting on local DoF '_key' at site indices [begin,end]
/// @param begin,end Iterator pair to list of site indices that define subset of
/// sites of interest
/// @param _syminfo SupercellSymInfo object that defines all symmetry properties
/// of supercell
/// @param _key DoFKey specifying which local DoF is of interest
/// @param _group vector of PermuteIterators forming the group that is to be
/// represented (this may be larger than a crystallographic factor group)
/// \result Orthogonal matrix comprising normal coordinate basis vectors as its
/// rows.
template <typename IterType>
Eigen::MatrixXd collective_dof_normal_coords(
    IterType begin, IterType end, SupercellSymInfo const &_syminfo,
    DoFKey const &_key, std::vector<PermuteIterator> const &_group);

/// \brief Find symmetry-adapted normal coordinate basis vectors for action of
/// '_group' acting on local DoF '_key' at site indices [begin,end] In addition
/// to returning the normal coordinate transformation matrix, it also returns a
/// list of dimensions of the corresponding irreducible invariant subspaces in
/// which the normal coordinate basis vectors reside
/// @param begin,end Iterator pair to list of site indices that define subset of
/// sites of interest
/// @param _syminfo SupercellSymInfo object that defines all symmetry properties
/// of supercell
/// @param _key DoFKey specifying which local DoF is of interest
/// @param _group vector of PermuteIterators forming the group that is to be
/// represented (this may be larger than a crystallographic factor group)
/// \result vector of IrrepInfo corresponding to irrep decomposition of the
/// action of the provided group on the specified DoF type, within the provided
/// supercell
template <typename IterType>
std::vector<SymRepTools::IrrepInfo> irrep_decomposition(
    IterType begin, IterType end, SupercellSymInfo const &_syminfo,
    DoFKey const &_key, std::vector<PermuteIterator> const &_group);

/// \brief Find symmetry-adapted normal coordinate basis vectors for action of
/// '_group' acting on local DoF '_key' at site indices [begin,end] In addition
/// to returning the normal coordinate transformation matrix, it also returns a
/// list of dimensions of the corresponding irreducible invariant subspaces in
/// which the normal coordinate basis vectors reside
/// @param begin,end Iterator pair to list of site indices that define subset of
/// sites of interest
/// @param _syminfo SupercellSymInfo object that defines all symmetry properties
/// of supercell
/// @param _key DoFKey specifying which local DoF is of interest
/// @param _group vector of PermuteIterators forming the group that is to be
/// represented (this may be larger than a crystallographic factor group)
/// @param __subspace matrix whose columns span a subspace of allowed DoF values
/// for selected sites \result vector of IrrepInfo corresponding to irrep
/// decomposition of the action of the provided group on the specified DoF type,
/// within the provided supercell
///         and restricted to the specified subspace
template <typename IterType>
std::vector<SymRepTools::IrrepInfo> irrep_decomposition(
    IterType begin, IterType end, SupercellSymInfo const &_syminfo,
    DoFKey const &_key, std::vector<PermuteIterator> const &_group,
    Eigen::Ref<const Eigen::MatrixXd> const &_subspace);

}  // namespace CASM

#endif
