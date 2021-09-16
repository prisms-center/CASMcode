#ifndef CASM_config_SupercellSymOp_impl
#define CASM_config_SupercellSymOp_impl

#include "casm/clexulator/ConfigDoFValues.hh"
#include "casm/configuration/PrimSymInfo.hh"
#include "casm/configuration/SupercellSymInfo.hh"
#include "casm/configuration/SupercellSymOp.hh"
#include "casm/crystallography/SymType.hh"

namespace CASM {
namespace config {

// --- Inline definitions ---

/// Default invalid SupercellSymOp, not equal to end iterator
SupercellSymOp::SupercellSymOp()
    : m_factor_group_index(), m_translation_index() {}

/// Construct SupercellSymOp
SupercellSymOp::SupercellSymOp(
    std::shared_ptr<Supercell const> const &_supercell,
    Index _factor_group_index, Index _translation_index)
    : m_supercell(_supercell),
      m_factor_group_index(_factor_group_index),
      m_translation_index(_translation_index),
      m_N_translation_index(
          m_supercell->sym_info.translation_permutations.size()) {}

/// Construct iterator
SupercellSymOp(std::shared_ptr<Supercell const> const &_supercell,
               Index _factor_group_index, Index _translation_index);

inline std::shared_ptr<Supercell const> const &SupercellSymOp::supercell()
    const {
  return m_supercell;
}

inline Index SupercellSymOp::factor_group_index() const {
  return m_factor_group_index;
}

inline Index SupercellSymOp::translation_index() const {
  return m_translation_index;
}

/// Returns a reference to this -- allows SupercellSymOp to be treated as an
/// iterator to SupercellSymOp object
inline SupercellSymOp const &SupercellSymOp::operator*() const { return *this; }

/// Returns a pointer to this -- allows SupercellSymOp to be treated as an
/// iterator to SupercellSymOp object
inline SupercellSymOp const *SupercellSymOp::operator->() const { return this; }

/// \brief prefix ++SupercellSymOp
inline SupercellSymOp &SupercellSymOp::operator++() {
  m_translation_index++;
  if (m_translation_index == m_N_translation) {
    m_translation_index = 0;
    m_factor_group_index++;
  }
  return *this;
}

/// \brief postfix SupercellSymOp++
inline SupercellSymOp SupercellSymOp::operator++(int) {
  SupercellSymOp cp(*this);
  ++cp;
  return cp;
}

/// \brief prefix --SupercellSymOp
inline SupercellSymOp &SupercellSymOp::operator--() {
  if (m_translation_index == 0) {
    m_factor_group_index--;
    m_translation_index = m_N_translation;
  }
  m_translation_index--;
  return *this;
}

/// \brief postfix SupercellSymOp--
inline SupercellSymOp SupercellSymOp::operator--(int) {
  PermuteIterator cp(*this);
  --cp;
  return cp;
}

/// \brief Return the SymOp for the current operation
///
/// Defined by:
///
///   translation_op * factor_group_op
///
/// In other words, the symmetry operation equivalent to application of the
/// factor group operation, FOLLOWED BY application of the translation
/// operation;
SymOp SupercellSymOp::sym_op() const {
  UnitCell translation_frac =
      this->m_supercell->unitcell_index_converter(this->m_translation_index);

  Eigen::Matrix3d const &prim_lat_column_mat =
      this->m_supercell->superlattice.prim_lattice().lat_column_mat();

  SymOp translation_cart = prim_lat_column_mat * translation_frac.cast<double>;

  SymOp const &fg_op =
      *this->m_supercell->factor_group.element[m_factor_group_index];

  return SymOp{fg_op.matrix, translation_cart + fg_op.translation,
               fg_op.is_time_reversal_active};
}

/// Returns the combination of factor group operation permutation and
/// translation permutation
Permutation SupercellSymOp::combined_permute() const {
  SupercellSymInfo const &sym_info = m_supercell.sym_info;
  return combined_permute(
      sym_info.factor_group_permutations[m_factor_group_index],  // first
      sym_info.translation_permutations[m_translation_index]);   // second
}

/// \brief Returns the inverse supercell operation
inline SupercellSymOp SupercellSymOp::inverse() const {
  // Copy *this, then update m_factor_group_index and m_translation_index
  SupercellSymOp inverse_it(*this);

  // Finding the inverse factor_group operation is straightforward
  Group const &supercell_factor_group = *this->m_supercell->factor_group;
  Index inverse_fg_index =
      supercell_factor_group.inverse_index[this->m_factor_group_index];
  inverse_it.m_factor_group_index = inverse_fg_index;

  // New translation can be found comparing the translation for the inverse of
  // the "total sym_op" of *this to the inverse of the untranslated sym_op

  // Find the translation (cartesian coordinates)
  SymOp const &inverse_sym_op = inverse(this->sym_op());
  SymOp const &inverse_fg_op = supercell_factor_group.element[inverse_fg_index];
  Eigen::Vector3d translation_cart =
      (inverse_sym_op.translation - inverse_fg_op.translation);

  // convert to fractional coordinates
  Superlattice const &superlattice = this->m_supercell->superlattice;
  UnitCell translation_uc =
      UnitCell::from_cartesian(translation_cart, superlattice.prim_lattice());

  // convert to linear index
  inverse_it.m_translation_index =
      m_supercell->unitcell_index_converter(translation_uc);

  return inverse_it;
}

/// \brief Returns the supercell operation equivalent to applying first RHS
/// and then *this
inline SupercellSymOp SupercellSymOp::operator*(
    SupercellSymOp const &RHS) const {
  // Copy *this, then update m_factor_group_index and m_translation_index
  PermuteIterator product_it(*this);

  // Finding the factor_group product is straightforward
  Group const &supercell_factor_group = *this->m_supercell->factor_group;
  product_it.m_factor_group_index =
      supercell_factor_group.multiplication_table[this->m_factor_group_index]
                                                 [RHS.m_factor_group_index];

  // New translation can be found comparing the translation for the product of
  // the "total sym_op" to the just the product factor group op translation

  // Find the translation (cartesian coordinates)
  SymOp total_product_op = this->sym_op() * RHS.sym_op();
  SymOp const &product_fg_op =
      supercell_factor_group
          .element[product_it.m_factor_group_index] Eigen::Vector3d
              translation_cart =
          total_product_op.translation - product_fg_op.translation;

  // convert to fractional coordinates
  Superlattice const &superlattice = this->m_supercell->superlattice;
  UnitCell translation_uc =
      UnitCell::from_cartesian(translation_cart, superlattice.prim_lattice());

  // convert to linear index
  product_it.m_translation_index =
      m_supercell->unitcell_index_converter(translation_uc);

  return product_it;
}

/// \brief Less than comparison (used to implement operator<() and other
/// standard comparisons via Comparisons)
inline bool SupercellSymOp::operator<(SupercellSymOp const &RHS) const {
  if (this->m_factor_group_index == RHS.m_factor_group_index) {
    return this->m_translation_index < RHS.m_translation_index;
  }
  return this->m_factor_group_index < RHS.m_factor_group_index;
}

/// \brief Equality comparison (used to implement operator==)
inline bool SupercellSymOp::eq_impl(SupercellSymOp const &RHS) const {
  if (m_supercell == RHS.m_supercell &&
      m_factor_group_index == RHS.m_factor_group_index &&
      m_translation_index == RHS.m_translation_index) {
    return true;
  }
  return false;
}

/// \brief Return inverse SymOp
inline SymOp inverse(SymOp const &op) {
  // x' = R * x + T
  // R.inv * x' = x * R.inv * T
  // R.inv * x' - R.inv * T = x

  // SymOp matrix is unitary, so inverse is equivalent to transpose.

  return SymOp{op.matrix.transpose(), -(op.matrix.transpose() * op.translation),
               op.is_time_reversal_active};
}

/// \brief Apply a symmetry operation specified by a SupercellSymOp to
/// ConfigDoFValues
inline ConfigDoFValues &apply(SupercellSymOp const &op,
                              ConfigDoFValues &dof_values) {
  Supercell const &supercell = *op.supercell();
  SupercellSymInfo const &supercell_sym_info = supercell.sym_info;
  Prim const &prim = op.supercell()->prim;
  PrimSymInfo const &prim_sym_info = prim.sym_info;
  Index n_vol = supercell.superlattice.size();
  Index n_sublat = prim.basicstructure.basis().size();
  Index n_sites = n_vol * n_sublat;

  Index supercell_fg_index = op.factor_group_index();
  Index prim_fg_index =
      supercell_sym_info.factor_group.head_group_index[supercell_fg_index];

  for (auto &dof : m_global_dofs) {
    Eigen::MatrixXd const &M =
        prim_sym_info.global_dof_symgroup_rep.at(dof.first)[prim_fg_index];
    dof.second = M * dof.second;
  }

  Permutation combined_permute{op.combined_permute()};

  if (dof_values.occupation.size()) {
    // permute occupant indices (if anisotropic)
    Eigen::VectorXi tmp{dof_values.occupation};
    if (prim_sym_info.has_aniso_occs) {
      Index l = 0;
      for (Index b = 0; b < n_sublat; ++b) {
        for (Index n = 0; n < n_vol; ++n, ++l) {
          Permutation const &occ_perm =
              prim_sym_info.occ_symgroup_rep[prim_fg_index][b];
          tmp[l] = occ_perm[tmp[l]];
        }
      }
    }

    // permute values amongst sites
    for (Index l = 0; l < n_sites; ++l) {
      dof_values.occupation[l] = tmp[combined_permute[l]];
    }
  }

  using clexulator::sublattice_block;
  for (auto &dof : m_dof_values.local_dof_values) {
    // vector of matrix, one per sublattice
    LocalDoFSymopRep const &local_dof_symop_rep =
        prim_sym_info.local_dof_symgroup_rep.at(dof.first)[prim_fg_index];

    // transform values on initial sites
    Eigen::MatrixXd const &init_value = dof.second;
    Eigen::MatrixXd tmp{init_value};
    for (Index b = 0; b < n_sublat; ++b) {
      Eigen::MatrixXd const &M = local_dof_symop_rep[b];
      Index dim = M.cols();
      sublattice_block(tmp, b, n_vol).topRows(dim) =
          M * sublattice_block(init_value, b, n_vol).topRows(dim);
    }

    // permute values amongst sites
    for (Index l = 0; l < n_sites; ++l) {
      dof.second.col(l) = tmp.col(combined_permute[l]);
    }
  }

  return *this;
}

/// \brief Apply a symmetry operation specified by a SupercellSymOp to
/// ConfigDoFValues
inline ConfigDoFValues copy_apply(SupercellSymOp const &op,
                                  ConfigDoFValues dof_values) {
  apply(op, dof_values);
  return dof_values;
}

}  // namespace config
}  // namespace CASM

#endif
