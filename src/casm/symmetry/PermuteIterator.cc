#include "casm/symmetry/PermuteIterator.hh"

#include <vector>

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/external/Eigen/src/Core/Matrix.h"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

PermuteIterator::PermuteIterator()
    : m_tmp_translation_permute(0), m_tmp_translation_index(-1) {}

PermuteIterator::PermuteIterator(const PermuteIterator &iter)
    : m_sym_info(iter.m_sym_info),
      m_factor_group_index(iter.m_factor_group_index),
      m_translation_index(iter.m_translation_index),
      m_tmp_translation_permute(iter.m_tmp_translation_permute),
      m_tmp_translation_index(iter.m_tmp_translation_index) {}

PermuteIterator::PermuteIterator(SupercellSymInfo const &_sym_info,
                                 Index _factor_group_index,
                                 Index _translation_index)
    : m_sym_info(&_sym_info),
      m_factor_group_index(_factor_group_index),
      m_translation_index(_translation_index),
      m_tmp_translation_permute(0),
      m_tmp_translation_index(-1) {}

PermuteIterator &PermuteIterator::operator=(PermuteIterator iter) {
  swap(*this, iter);
  return *this;
}

/// Returns a reference to this
const PermuteIterator &PermuteIterator::operator*() const { return *this; }

/// Returns a pointer to this
const PermuteIterator *PermuteIterator::operator->() const { return this; }

/// Returns the combination of factor_group permutation and translation
/// permutation
Permutation PermuteIterator::combined_permute() const {
  return translation_permute() * factor_group_permute();
}

/// Reference the SupercellSymInfo containing the operations being pointed at
SupercellSymInfo const &PermuteIterator::sym_info() const {
  return *m_sym_info;
}

/// Reference the supercell factor group (`SupercellSymInfo::factor_group()`)
SymGroup const &PermuteIterator::factor_group() const {
  return sym_info().factor_group();
}

/// Return the supercell factor group index
///
/// Index into `SupercellSymInfo::factor_group()` of the factor group op being
/// pointed at
Index PermuteIterator::factor_group_index() const {
  return m_factor_group_index;
}

/// Return the prim factor group index
///
/// Index into `Structure::factor_group()` of the factor group op being pointed
/// at
Index PermuteIterator::prim_factor_group_index() const {
  return factor_group()[factor_group_index()].master_group_index();
}

/// Return the index into the supercell translation permutations
///
/// Index into `SupercellSymInfo::translation_permutations()` of the translation
/// being pointed at
Index PermuteIterator::translation_index() const { return m_translation_index; }

/// Return the factor group permutation being pointed at
const Permutation &PermuteIterator::factor_group_permute() const {
  return *(sym_info()
               .site_permutation_symrep()[m_factor_group_index]
               ->permutation());
}

/// Return the translation permutation being pointed at
const Permutation &PermuteIterator::translation_permute() const {
  if (m_sym_info->translation_permutations().size() != 0) {
    return m_sym_info->translation_permutations()[m_translation_index];
  }
  if (m_translation_index != m_tmp_translation_index) {
    m_tmp_translation_index = m_translation_index;
    m_tmp_translation_permute = make_translation_permutation(
        m_tmp_translation_index, m_sym_info->unitcellcoord_index_converter(),
        m_sym_info->unitcell_index_converter());
  }
  return m_tmp_translation_permute;
}

/// Returns representation of current operation corresponding to species
/// transformation on sublattice b
SymOpRepresentation const &PermuteIterator::occ_rep(Index b) const {
  return *sym_info().occ_symreps()[b][factor_group_index()];
}

/// Returns representation of current operation
/// corresponding to local DoF specified by _key on sublattice b
SymOpRepresentation const &PermuteIterator::local_dof_rep(DoFKey const &_key,
                                                          Index b) const {
  return *sym_info().local_dof_symreps(_key)[b][factor_group_index()];
}

/// Check if local dof rep is empty
bool PermuteIterator::local_dof_rep_empty(DoFKey const &_key, Index b) const {
  if (sym_info().local_dof_symreps(_key)[b].rep_ptr() == nullptr) {
    return true;
  }
  return false;
}

/// Returns representation of current operation corresponding to global DoF
/// specified by _key
SymOpRepresentation const &PermuteIterator::global_dof_rep(
    DoFKey const &_key) const {
  return *sym_info().global_dof_symrep(_key)[factor_group_index()];
}

SymOp PermuteIterator::sym_op() const {
  UnitCell translation_lattice_site =
      this->sym_info().unitcell_index_converter()(this->translation_index());
  SymOp lattice_translation_op = SymOp::translation(
      make_superlattice_coordinate(translation_lattice_site,
                                   this->sym_info().superlattice())
          .cart());
  return lattice_translation_op *
         sym_info().factor_group()[m_factor_group_index];
}

Index PermuteIterator::permute_ind(Index i) const {
  return factor_group_permute()[translation_permute()[i]];
}

bool PermuteIterator::operator<(const PermuteIterator &iter) const {
  if (this->factor_group_index() == iter.factor_group_index()) {
    return this->translation_index() < iter.translation_index();
  }
  return this->factor_group_index() < iter.factor_group_index();
}

bool PermuteIterator::eq_impl(const PermuteIterator &iter) const {
  if (m_sym_info == iter.m_sym_info &&
      m_factor_group_index == iter.m_factor_group_index &&
      m_translation_index == iter.m_translation_index)
    return true;
  return false;
}

// prefix ++PermuteIterator
PermuteIterator &PermuteIterator::operator++() {
  m_translation_index++;
  if (m_translation_index == m_sym_info->superlattice().size()) {
    m_translation_index = 0;
    m_factor_group_index++;
  }
  return *this;
}

// postfix PermuteIterator++
PermuteIterator PermuteIterator::operator++(int) {
  PermuteIterator cp(*this);
  ++cp;
  return cp;
}

// prefix --PermuteIterator
PermuteIterator &PermuteIterator::operator--() {
  if (m_translation_index == 0) {
    m_factor_group_index--;
    m_translation_index = m_sym_info->superlattice().size();
  }
  m_translation_index--;
  return *this;
}

// postfix PermuteIterator--
PermuteIterator PermuteIterator::operator--(int) {
  PermuteIterator cp(*this);
  --cp;
  return cp;
}

PermuteIterator PermuteIterator::begin_next_fg_op() const {
  PermuteIterator it(*this);
  it.m_translation_index = 0;
  it.m_factor_group_index++;
  return it;
}

// Might be able to do this in a faster way, but would be much harder to
// understand
PermuteIterator PermuteIterator::inverse() const {
  PermuteIterator it(*this);
  // Finding the inverse factor_group operation is straightforward
  it.m_factor_group_index = factor_group().ind_inverse(factor_group_index());

  // Easiest way to get the new translation is just to compare the tau of the
  // inverse of the 'total' sym_op (described by *this), to the inverse of the
  // untranslated symop (described by
  // m_fg_site_permutation_symrep.sym_op(it.m_factor_group_index)) Result is the
  // portion of the inverse sym_op that needs to be described by a lattice point
  // translation
  Eigen::Vector3d translation_cart =
      (sym_op().inverse().tau() -
       factor_group()[it.factor_group_index()].tau());
  UnitCell translation_uc = UnitCell::from_cartesian(
      translation_cart, this->sym_info().prim_lattice());
  it.m_translation_index =
      this->sym_info().unitcell_index_converter()(translation_uc);

  return it;
}

PermuteIterator PermuteIterator::operator*(const PermuteIterator &RHS) const {
  PermuteIterator it(*this);
  // Finding the inverse factor_group operation is straightforward
  it.m_factor_group_index =
      factor_group().ind_prod(factor_group_index(), RHS.factor_group_index());

  // Easiest way to get the new translation is just to compare the tau of the
  // 'total' sym_op (described by (*this).sym_op()*RHS.sym_op()), to the
  // untranslated symop product (described by
  // m_fg_site_permutation_symrep.sym_op(it.factor_group_index())) Result is the
  // portion of the product sym_op that needs to be described by a lattice point
  // translation
  Eigen::Vector3d translation_cart =
      (sym_op() * RHS.sym_op()).tau() -
      factor_group()[it.factor_group_index()].tau();
  UnitCell translation_uc = UnitCell::from_cartesian(
      translation_cart, this->sym_info().prim_lattice());
  it.m_translation_index =
      this->sym_info().unitcell_index_converter()(translation_uc);

  return it;
}

/// \brief Returns a SymGroup generated from a range of PermuteIterator
///
/// \param begin,end A range of PermuteIterator
///
/// - The result is sorted
/// - The result uses the Supercell lattice for periodic comparisons
template <typename PermuteIteratorIt>
SymGroup make_point_group(PermuteIteratorIt begin, PermuteIteratorIt end,
                          const Lattice &supercell_lattice) {
  SymGroup result;
  result.set_lattice(supercell_lattice);
  for (; begin != end; ++begin) {
    Index f = begin->sym_op().index();
    Index i;
    for (i = 0; i < result.size(); ++i) {
      if (f == result[i].index()) break;
    }
    if (i == result.size()) {
      result.push_back((begin->sym_op()).no_trans());
      // std::cout << "pushed back op " << result.back().index() << "\n";
    }
  }
  result.sort();
  return result;
}

/// \brief Returns a SymGroup generated from a range of PermuteIterator
///
/// \param begin,end A range of PermuteIterator
///
/// - The result is sorted
/// - The result uses the Supercell lattice for periodic comparisons
template <typename PermuteIteratorIt>
SymGroup make_sym_group(PermuteIteratorIt begin, PermuteIteratorIt end,
                        const Lattice &supercell_lattice) {
  SymGroup sym_group;
  sym_group.set_lattice(supercell_lattice);
  while (begin != end) {
    sym_group.push_back(begin->sym_op());
    ++begin;
  }
  sym_group.sort();
  return sym_group;
}

/// \brief Returns a std::unique_ptr<SymGroup> generated from a range of
/// PermuteIterator
///
/// \param begin,end A range of PermuteIterator
///
/// - The result is sorted
/// - The result uses the Supercell lattice for periodic comparisons
template <typename PermuteIteratorIt>
std::unique_ptr<SymGroup> make_unique_sym_group(
    PermuteIteratorIt begin, PermuteIteratorIt end,
    const Lattice &supercell_lattice) {
  auto sym_group_ptr = notstd::make_unique<SymGroup>();
  sym_group_ptr->set_lattice(supercell_lattice);
  while (begin != end) {
    sym_group_ptr->push_back(begin->sym_op());
    ++begin;
  }
  sym_group_ptr->sort();
  return sym_group_ptr;
}

template SymGroup make_sym_group(PermuteIterator begin, PermuteIterator end,
                                 const Lattice &supercell_lattice);
template SymGroup make_sym_group(
    std::vector<PermuteIterator>::const_iterator begin,
    std::vector<PermuteIterator>::const_iterator end,
    const Lattice &supercell_lattice);
template SymGroup make_sym_group(std::vector<PermuteIterator>::iterator begin,
                                 std::vector<PermuteIterator>::iterator end,
                                 const Lattice &supercell_lattice);

template SymGroup make_point_group(PermuteIterator begin, PermuteIterator end,
                                   const Lattice &supercell_lattice);
template SymGroup make_point_group(
    std::vector<PermuteIterator>::const_iterator begin,
    std::vector<PermuteIterator>::const_iterator end,
    const Lattice &supercell_lattice);
template SymGroup make_point_group(std::vector<PermuteIterator>::iterator begin,
                                   std::vector<PermuteIterator>::iterator end,
                                   const Lattice &supercell_lattice);

/// \brief Make permute group for local property symmetry in a supercell
///
/// \brief generating_group Local property group. Must be a local property
///     group, in which each prim factor group operation only appears
/// \brief factor_group The prim factor group.
/// \brief supercell_sym_info SupercellSymInfo for the supercell in which
///     the local permute group is being generated.
///
/// \returns local_permute_group, the PermuteIterator consistent with both
///     the supercell and the local point group
std::vector<PermuteIterator> make_local_permute_group(
    const SymGroup &local_generating_group, const SymGroup &factor_group,
    SupercellSymInfo const &supercell_sym_info) {
  std::set<Index> generating_group_indices;
  for (SymOp const &op : local_generating_group) {
    generating_group_indices.emplace(op.index());
  }

  std::map<Index, Index> prim_to_super_factor_group_index;
  Index supercell_factor_group_index = 0;
  for (SymOp const &op : supercell_sym_info.factor_group()) {
    prim_to_super_factor_group_index[op.index()] = supercell_factor_group_index;
    ++supercell_factor_group_index;
  }

  std::vector<PermuteIterator> result;
  try {
    for (SymOp const &op : local_generating_group) {
      if (prim_to_super_factor_group_index.count(op.index())) {
        Eigen::Vector3d trans = op.tau() - factor_group[op.index()].tau();
        PermuteIterator permute = supercell_sym_info.permute_it(
            prim_to_super_factor_group_index.at(op.index()),
            UnitCell::from_cartesian(trans, supercell_sym_info.prim_lattice()));
        result.push_back(permute);
      }
    }
  } catch (std::exception &e) {
    std::stringstream msg;
    msg << "Error in permute_group: failed to generate PermuteIterator: "
        << e.what();
    throw e;
  }

  return result;
}

/// \brief Filter PermuteIterator to keep only operations consistent
///     with a factor group in a sub-supercell
///
/// \param supercell_permute_begin,supercell_permute_end A range of
///     PermuteIterator describing the symmetry in a supercell
/// \param subsupercell_factor_group A set of PermuteIterator in a
///     sub-supercell, specifying a factor group in a sub-supercell.
///
/// \returns allowed_permute, a vector of PermuteIterator containing
///     the operations in [supercell_permute_begin,supercell_permute_end]
///     that are equivalent to one of the operations in subsupercell_permute
///     up to translations of the sub-supercell lattice.
template <typename PermuteIteratorIt>
std::vector<PermuteIterator> make_allowed_permute(
    PermuteIteratorIt supercell_permute_begin,
    PermuteIteratorIt supercell_permute_end,
    std::set<PermuteIterator> const &subsupercell_factor_group) {
  std::vector<PermuteIterator> result;
  if (!subsupercell_factor_group.size()) {
    return result;
  }
  if (supercell_permute_begin == supercell_permute_end) {
    return result;
  }
  auto const &subsupercell_sym_info =
      subsupercell_factor_group.begin()->sym_info();
  auto const &supercell_sym_info = supercell_permute_begin->sym_info();

  if (!is_superlattice(supercell_sym_info.supercell_lattice(),
                       subsupercell_sym_info.supercell_lattice(),
                       subsupercell_sym_info.supercell_lattice().tol())
           .first) {
    return result;
  }
  auto const &subsupercell_converter =
      subsupercell_sym_info.unitcell_index_converter();
  auto const &supercell_converter =
      supercell_sym_info.unitcell_index_converter();
  for (auto it = supercell_permute_begin; it != supercell_permute_end; ++it) {
    PermuteIterator test_permute = subsupercell_sym_info.permute_it(
        it->factor_group_index(),
        subsupercell_converter(supercell_converter(it->translation_index())));
    if (subsupercell_factor_group.count(test_permute)) {
      result.push_back(*it);
    }
  }
  return result;
}

template std::vector<PermuteIterator> make_allowed_permute(
    PermuteIterator supercell_permute_begin,
    PermuteIterator supercell_permute_end,
    std::set<PermuteIterator> const &subsupercell_permute);

template std::vector<PermuteIterator> make_allowed_permute(
    std::vector<PermuteIterator>::const_iterator supercell_permute_begin,
    std::vector<PermuteIterator>::const_iterator supercell_permute_end,
    std::set<PermuteIterator> const &subsupercell_permute);

template std::vector<PermuteIterator> make_allowed_permute(
    std::vector<PermuteIterator>::iterator supercell_permute_begin,
    std::vector<PermuteIterator>::iterator supercell_permute_end,
    std::set<PermuteIterator> const &subsupercell_permute);

void swap(PermuteIterator &a, PermuteIterator &b) {
  std::swap(a.m_sym_info, b.m_sym_info);
  std::swap(a.m_factor_group_index, b.m_factor_group_index);
  std::swap(a.m_translation_index, b.m_translation_index);
  std::swap(a.m_tmp_translation_permute, b.m_tmp_translation_permute);
  std::swap(a.m_tmp_translation_index, b.m_tmp_translation_index);
}

/// Return true if the permutation does not given sites and other sites
bool site_indices_are_invariant(PermuteIterator const &permute_it,
                                std::set<Index> const &site_indices) {
  // Applying the permutation indicated by `permute_it` moves the value from
  // site index `permute_it.permute_ind(s)` to site index `s`, for each `s` in
  // the set. Therefore, if none of `permute_it.permute_ind(s)` are outside the
  // set `site_indices` the sites are invariant.

  return std::none_of(site_indices.begin(), site_indices.end(), [&](Index s) {
    return site_indices.count(permute_it.permute_ind(s)) == 0;
  });
}

jsonParser &to_json(const PermuteIterator &it, jsonParser &json) {
  json.put_obj();
  json["factgrp"] = it.factor_group_index();
  json["trans"] = it.translation_index();
  return json;
}

PermuteIterator jsonConstructor<PermuteIterator>::from_json(
    const jsonParser &json, const SupercellSymInfo &scel_info) {
  return scel_info.permute_it(json["factgrp"].get<Index>(),
                              json["trans"].get<Index>());
}

namespace adapter {

SymOp Adapter<SymOp, PermuteIterator>::operator()(
    PermuteIterator const &adaptable) const {
  return adaptable.sym_op();
}
}  // namespace adapter

}  // namespace CASM
