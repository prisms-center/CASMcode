#include "casm/symmetry/SupercellSymInfo.hh"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include "casm/casm_io/container/stream_io.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/LatticeIsEquivalent.hh"
#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/crystallography/SymTools.hh"
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

/// \brief Make a single translation permutation
///
/// \param translation_index Index in range [0, n_unitcell)
/// \param bijk_index_converter Site index converter
/// \param ijk_index_converter Unitcell index converter
///
/// \returns permutation, A permutation that translates all
///     supercell occupants by a multiple of the lattice vectors
Permutation make_translation_permutation(
    Index translation_index,
    xtal::UnitCellCoordIndexConverter bijk_index_converter,
    xtal::UnitCellIndexConverter ijk_index_converter) {
  std::vector<Index> single_translation_permutation(
      bijk_index_converter.total_sites(), -1);
  UnitCell translation_uc = ijk_index_converter(translation_index);

  // Loops over all the sites
  for (Index old_site_ix = 0; old_site_ix < bijk_index_converter.total_sites();
       ++old_site_ix) {
    UnitCellCoord old_site_ucc = bijk_index_converter(old_site_ix);
    Index new_site_ix = bijk_index_converter(old_site_ucc + translation_uc);

    single_translation_permutation[new_site_ix] = old_site_ix;
  }
  // You should have given a permutation value to every single site
  assert(std::find(single_translation_permutation.begin(),
                   single_translation_permutation.end(),
                   -1) == single_translation_permutation.end());
  return Permutation(single_translation_permutation);
}

/// \brief Make all translation permutations
std::vector<Permutation> make_translation_permutations(
    const Eigen::Matrix3l &transformation_matrix, int basis_sites_in_prim) {
  xtal::UnitCellCoordIndexConverter bijk_index_converter(transformation_matrix,
                                                         basis_sites_in_prim);
  xtal::UnitCellIndexConverter ijk_index_converter(transformation_matrix);
  std::vector<Permutation> translation_permutations;

  // Loops over lattice points
  for (Index translation_ix = 0;
       translation_ix < ijk_index_converter.total_sites(); ++translation_ix) {
    translation_permutations.push_back(make_translation_permutation(
        translation_ix, bijk_index_converter, ijk_index_converter));
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
      m_factor_group(sym::invariant_subgroup(_prim_factor_group, _super_lat)),
      m_basis_perm_symrep(factor_group(), basis_permutation_symrep_ID),
      m_has_aniso_occs(false),
      m_has_occupation_dofs(false) {
  if (m_supercell_superlattice.size() <= 100) {
    m_translation_permutations = make_translation_permutations(
        this->superlattice().transformation_matrix_to_super(),
        num_sites_in_prim);
  }

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

std::string hermite_normal_form_name(const Eigen::Matrix3l &matrix) {
  std::string name_str;

  Eigen::Matrix3i H = hermite_normal_form(matrix.cast<int>()).first;
  name_str = "SCEL";
  std::stringstream tname;
  // Consider using a for loop with HermiteCounter_impl::_canonical_unroll here
  tname << H(0, 0) * H(1, 1) * H(2, 2) << "_" << H(0, 0) << "_" << H(1, 1)
        << "_" << H(2, 2) << "_" << H(1, 2) << "_" << H(0, 2) << "_" << H(0, 1);
  name_str.append(tname.str());

  return name_str;
}

Eigen::Matrix3l make_hermite_normal_form(std::string hermite_normal_form_name) {
  std::vector<std::string> tmp, tokens;
  try {
    // else construct transf_mat from name (make sure to remove any empty
    // tokens)
    boost::split(tmp, hermite_normal_form_name, boost::is_any_of("SCEL_"),
                 boost::token_compress_on);
    std::copy_if(tmp.begin(), tmp.end(), std::back_inserter(tokens),
                 [](const std::string &val) { return !val.empty(); });
    if (tokens.size() != 7) {
      throw std::invalid_argument(
          "Error in make_supercell: supercell name format error");
    }
    Eigen::Matrix3l T;
    auto cast = [](std::string val) { return std::stol(val); };
    T << cast(tokens[1]), cast(tokens[6]), cast(tokens[5]), 0, cast(tokens[2]),
        cast(tokens[4]), 0, 0, cast(tokens[3]);
    return T;
  } catch (std::exception &e) {
    std::string format = "SCELV_T00_T11_T22_T12_T02_T01";
    std::stringstream ss;
    ss << "Error in make_hermite_normal_form: "
       << "expected format: " << format << ", "
       << "name: |" << hermite_normal_form_name << "|"
       << ", "
       << "tokens: " << tokens << ", "
       << "tokens.size(): " << tokens.size() << ", "
       << "error: " << e.what();
    throw std::runtime_error(ss.str());
  }
}

/// Make the supercell name from SupercellSymInfo
///
/// For supercells that are equivalent to the canonical supercell:
/// - The supercell name is `SCELV_A_B_C_D_E_F`, where 'V' is supercell volume
///   (number of unit cells), and 'A-F' are the six non-zero elements of the
///   hermite normal form of the canonical supercell transformation matrix
///   (T00, T11, T22, T12, T02, T01)
///
/// For supercells that are not equivalent to the canonical supercell:
/// - The supercell name is `$CANON_SCELNAME.$FG_INDEX` where CANON_SCELNAME is
///   the supercell name for the canonical equivalent supercell and FG_INDEX is
///   the lowest index prim factor group operation which can be applied to the
///   canonical supercell lattice to construct a lattice equivalent to the
///   input supercell lattice.
///
std::string make_supercell_name(SymGroup const &point_group,
                                Lattice const &prim_lattice,
                                Lattice const &supercell_lattice) {
  xtal::Superlattice canon_superlattice{
      prim_lattice,
      xtal::canonical::equivalent(supercell_lattice, point_group)};
  std::string supercell_name = hermite_normal_form_name(
      canon_superlattice.transformation_matrix_to_super());
  if (!xtal::is_equivalent(supercell_lattice,
                           canon_superlattice.superlattice())) {
    double tol = prim_lattice.tol();
    auto from_canonical_index =
        is_equivalent_superlattice(supercell_lattice,
                                   canon_superlattice.superlattice(),
                                   point_group.begin(), point_group.end(), tol)
            .first->index();
    supercell_name += ("." + std::to_string(from_canonical_index));
  }
  return supercell_name;
}

/// Make the canonical supercell name from SupercellSymInfo
std::string make_canonical_supercell_name(SymGroup const &point_group,
                                          Lattice const &prim_lattice,
                                          Lattice const &supercell_lattice) {
  xtal::Superlattice canon_superlattice{
      prim_lattice,
      xtal::canonical::equivalent(supercell_lattice, point_group)};
  return hermite_normal_form_name(
      canon_superlattice.transformation_matrix_to_super());
}

/// Construct a Superlattice from the supercell name
xtal::Superlattice make_superlattice_from_supercell_name(
    SymGroup const &factor_group, Lattice const &prim_lattice,
    std::string supercell_name) {
  try {
    // tokenize name: check if non-canonical
    std::vector<std::string> tokens;
    boost::split(tokens, supercell_name, boost::is_any_of("."),
                 boost::token_compress_on);

    // validate name
    if (tokens.size() == 0 || tokens.size() > 2) {
      throw std::invalid_argument("supercell_name format error");
    }

    Eigen::Matrix3l T = make_hermite_normal_form(tokens[0]);
    xtal::Lattice super_lattice = make_superlattice(prim_lattice, T);

    if (tokens.size() == 2) {
      Index fg_op_index = std::stol(tokens[1]);
      if (fg_op_index >= factor_group.size()) {
        std::stringstream ss;
        ss << "Error in make_superlattice_from_supercell_name: "
           << "found prim factor group index: " << fg_op_index
           << ", which is out of range [0, " << factor_group.size() << ").";
        throw std::invalid_argument(ss.str());
      }
      super_lattice = sym::copy_apply(factor_group[fg_op_index], super_lattice);
      // ** uses super_lattice point group **
      super_lattice = xtal::canonical::equivalent(super_lattice);
    } else {
      // ** uses point group of provided factor group (typically from prim) **
      SymGroup point_group = factor_group.copy_no_trans();
      super_lattice = xtal::canonical::equivalent(super_lattice, point_group);
    }
    return xtal::Superlattice{prim_lattice, super_lattice};

  } catch (std::exception &e) {
    std::string format = "$CANON_SCEL_NAME[.$PRIM_FG_OP]";
    std::stringstream ss;
    ss << "Error in make_superlattice_from_supercell_name: "
       << "expected format: " << format << ", name: " << supercell_name
       << ", error: " << e.what();
    throw std::runtime_error(ss.str());
  }
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
