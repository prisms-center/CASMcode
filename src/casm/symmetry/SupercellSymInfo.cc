#include "casm/symmetry/SupercellSymInfo.hh"

#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/symmetry/SupercellSymInfo.hh"
#include "casm/symmetry/SymBasisPermute.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymMatrixXd.hh"
#include "casm/symmetry/SymPermutation.hh"
#include "casm/symmetry/SymTools.hh"

namespace CASM {

std::vector<Permutation> make_translation_permutations(
    const Eigen::Matrix3l &transformation_matrix, int basis_sites_in_prim) {
  xtal::UnitCellCoordIndexConverter bijk_index_converter(transformation_matrix,
                                                         basis_sites_in_prim);
  xtal::UnitCellIndexConverter ijk_index_converter(transformation_matrix);
  std::vector<Permutation> translation_permutations;

  // Loops over lattice points
  for (Index translation_ix = 0;
       translation_ix < ijk_index_converter.total_sites(); ++translation_ix) {
    std::vector<Index> single_translation_permutation(
        bijk_index_converter.total_sites(), -1);
    UnitCell translation_uc = ijk_index_converter(translation_ix);

    // Loops over all the sites
    for (Index old_site_ix = 0;
         old_site_ix < bijk_index_converter.total_sites(); ++old_site_ix) {
      UnitCellCoord old_site_ucc = bijk_index_converter(old_site_ix);
      Index new_site_ix = bijk_index_converter(old_site_ucc + translation_uc);

      single_translation_permutation[new_site_ix] = old_site_ix;
    }
    // You should have given a permutation value to every single site
    assert(std::find(single_translation_permutation.begin(),
                     single_translation_permutation.end(),
                     -1) == single_translation_permutation.end());
    translation_permutations.push_back(
        Permutation(single_translation_permutation));
  }
  return translation_permutations;
}

SupercellSymInfo::SupercellSymInfo(
    Lattice const &_prim_lat, Lattice const &_super_lat,
    Index num_sites_in_prim, SymGroup const &_prim_factor_group,
    SymGroupRepID basis_permutation_symrep_ID,
    std::map<DoFKey, SymGroupRepID> const &global_dof_symrep_IDs,
    std::vector<SymGroupRepID> const &occ_symrep_IDs,
    std::map<DoFKey, std::vector<SymGroupRepID> > const &local_dof_symrep_IDs)
    : m_supercell_superlattice(_prim_lat, _super_lat),
      m_unitcell_to_index_converter(
          m_supercell_superlattice.transformation_matrix_to_super()),
      m_unitcellcoord_to_index_converter(
          m_supercell_superlattice.transformation_matrix_to_super(),
          num_sites_in_prim),
      m_translation_permutations(make_translation_permutations(
          this->superlattice().transformation_matrix_to_super(),
          num_sites_in_prim)),
      m_factor_group(sym::invariant_subgroup(_prim_factor_group, _super_lat)),
      m_basis_perm_symrep(factor_group(), basis_permutation_symrep_ID),
      m_has_aniso_occs(false),
      m_has_occupation_dofs(false) {
  for (auto const &dofID : global_dof_symrep_IDs)
    m_global_dof_symreps.emplace(std::make_pair(
        dofID.first, SymGroupRep::RemoteHandle(factor_group(), dofID.second)));

  for (auto const &dofID : local_dof_symrep_IDs) {
    SublatSymReps treps(num_sites_in_prim);
    for (Index b = 0; b < num_sites_in_prim; ++b) {
      if (!dofID.second[b].empty())
        treps[b] = SymGroupRep::RemoteHandle(factor_group(), dofID.second[b]);
    }
    m_local_dof_symreps.emplace(std::make_pair(dofID.first, std::move(treps)));
  }

  m_occ_symreps.resize(num_sites_in_prim);
  for (Index b = 0; b < num_sites_in_prim; ++b) {
    if (!occ_symrep_IDs[b].is_identity()) {
      m_has_aniso_occs = true;
      m_has_occupation_dofs = true;
    } else if (occ_symrep_IDs[b].rep_index() > 1) {
      m_has_occupation_dofs = true;
    }

    m_occ_symreps[b] =
        SymGroupRep::RemoteHandle(factor_group(), occ_symrep_IDs[b]);
  }
}

SymGroupRepID make_permutation_representation(
    const SymGroup &group,
    const xtal::UnitCellCoordIndexConverter &bijk_index_converter,
    const Lattice &prim_lattice, const SymGroupRepID &prim_symrep_ID) {
  SymGroupRepID perm_rep_ID = group.allocate_representation();
  long total_sites = bijk_index_converter.total_sites();
  for (Index operation_ix = 0; operation_ix < group.size(); ++operation_ix) {
    const auto &operation = group[operation_ix];
    std::vector<Index> permutation(total_sites);
    for (Index old_l = 0; old_l < total_sites; ++old_l) {
      const UnitCellCoord &old_ucc = bijk_index_converter(old_l);
      UnitCellCoord new_ucc =
          sym::copy_apply(operation, old_ucc, prim_lattice, prim_symrep_ID);
      Index new_l = bijk_index_converter(new_ucc);
      permutation[new_l] = old_l;
    }

    group[operation_ix].set_rep(perm_rep_ID, SymPermutation(permutation));
  }
  return perm_rep_ID;
}

/// Returns a "RemoteHandle" to the sublattice permutation representation of the
/// supercell's factor group
///
/// A sublattice permutation consists of an index permutation and an integer
/// lattice translation
SymGroupRep::RemoteHandle const &SupercellSymInfo::basis_permutation_symrep()
    const {
  return m_basis_perm_symrep;
}

/// Returns a "RemoteHandle" to the site permutation representation
///
/// A site permutation consists of an index permutation for the sites of the
/// supercell. If
///
///     site_coordinate(index_after) = apply(operation,
///     site_coordinate(index_before)),
///
/// then this encodes permutations as
///
///     permutation[index_after] = index_before.
///
/// This convention treats fixed positions in space as having an unchanging
/// order, and the permutation describes a rearrangement of objects among those
/// sites, due to a transformation.
///
/// Note:
/// - The permutation representation is constructed for
/// `this->prim().factor_group()`, which may
///   contain more operations than `this->factor_group()`, so the Permutation
///   SymGroupRep may have 'gaps' at the operations that aren't in
///   `this->factor_group()`.
/// - This function returns the "RemoteHandle" to the permutation representation
/// which allows
///   accessing permutation representations using the
///   "supercell_factor_group_index" (index into `this->factor_group()`) instead
///   of the "prim_factor_group_index" (index into
///   `this->prim().factor_group()`).
///
SymGroupRep::RemoteHandle const &SupercellSymInfo::site_permutation_symrep()
    const {
  if (m_site_perm_symrep.empty()) {
    m_site_perm_symrep = SymGroupRep::RemoteHandle(
        this->factor_group(),
        make_permutation_representation(
            this->factor_group(), this->unitcellcoord_index_converter(),
            this->prim_lattice(),
            this->basis_permutation_symrep().symrep_ID()));
  }

  return m_site_perm_symrep;
}

/// Site permutation corresponding to supercell factor group operation
const Permutation &SupercellSymInfo::factor_group_permute(
    Index supercell_factor_group_index) const {
  return *(
      site_permutation_symrep()[supercell_factor_group_index]->permutation());
}

/// \brief Begin iterator over translation permutations
SupercellSymInfo::permute_const_iterator SupercellSymInfo::translate_begin()
    const {
  return permute_begin();
}

/// \brief End iterator over translation permutations
SupercellSymInfo::permute_const_iterator SupercellSymInfo::translate_end()
    const {
  return permute_begin().begin_next_fg_op();
}

SupercellSymInfo::permute_const_iterator SupercellSymInfo::permute_begin()
    const {
  return permute_it(0, 0);  // starting indices
}

SupercellSymInfo::permute_const_iterator SupercellSymInfo::permute_end() const {
  return permute_it(factor_group().size(), 0);  // one past final indices
}

SupercellSymInfo::permute_const_iterator SupercellSymInfo::permute_it(
    Index fg_index, Index trans_index) const {
  return permute_const_iterator(*this, fg_index, trans_index);
}

SupercellSymInfo::permute_const_iterator SupercellSymInfo::permute_it(
    Index fg_index, UnitCell trans) const {
  Index trans_index = this->unitcell_index_converter()(trans);
  return permute_it(fg_index, trans_index);
}

/// \brief Make the matrix representation for group '_group' describing the
/// transformation of DoF '_key' among a subset of sites
///
/// \param site_indices Set of site indices that define subset of sites of
///     interest
/// \param _syminfo SupercellSymInfo object that defines all symmetry properties
///     of supercell
/// \param _key DoFKey specifying which local DoF is of interest
/// \param _group vector of PermuteIterators forming the group that is to be
///     represented (this may be larger than a crystallographic factor group)
///
/// \result A std::pair containing a MasterSymGroup instantiation of _group and
/// a SymGroupRepID that can be used to access the 'collective_dof_symrep'
/// within the returned MasterSymGroup
///
std::pair<MasterSymGroup, SymGroupRepID> make_collective_dof_symrep(
    std::set<Index> const &site_indices, SupercellSymInfo const &_syminfo,
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
  for (Index site_index : site_indices) {
    Index b = _syminfo.unitcellcoord_index_converter()(site_index).sublattice();
    Index site_dof_dim = subreps[b].dim();
    site_index_to_basis_index[site_index] = total_dim;
    total_dim += site_dof_dim;
  }

  // make matrix rep, by filling in blocks with site dof symreps
  Eigen::MatrixXd trep(total_dim, total_dim);
  Index g = 0;
  for (PermuteIterator const &perm : _group) {
    trep.setZero();
    for (Index site_index : site_indices) {
      // "to_site" (after applying symmetry) determines row of block
      // can't fail, because it was built from [begin, end)
      Index to_site_index = site_index;
      Index row = site_index_to_basis_index.find(to_site_index)->second;

      // "from_site" (before applying symmetry) determines col of block
      // could fail, if mismatch between [begin, end) and group
      Index from_site_index = perm.permute_ind(site_index);
      auto col_it = site_index_to_basis_index.find(from_site_index);
      if (col_it == site_index_to_basis_index.end()) {
        throw std::runtime_error(
            "Error in collective_dof_symrep: Input group includes permutations "
            "between selected and unselected sites.");
      }
      Index col = col_it->second;

      // "from_site" sublattice and factor group op index
      // are used to lookup the site dof rep matrix
      Index from_site_b =
          _syminfo.unitcellcoord_index_converter()(from_site_index)
              .sublattice();
      Eigen::MatrixXd U =
          *(subreps[from_site_b][perm.factor_group_index()]->MatrixXd());

      // insert matrix as block in collective dof symrep
      trep.block(row, col, U.rows(), U.cols()) = U;
    }
    result.first[g++].set_rep(result.second, SymMatrixXd(trep));
  }
  result.first.sort();
  return result;
}

}  // namespace CASM
