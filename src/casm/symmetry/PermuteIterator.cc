#include "casm/symmetry/PermuteIterator.hh"

#include <vector>

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/external/Eigen/src/Core/Matrix.h"

namespace CASM {

PermuteIterator::PermuteIterator() {}

PermuteIterator::PermuteIterator(const PermuteIterator &iter)
    : m_sym_info(iter.m_sym_info),
      m_trans_permute(&(m_sym_info->translation_permutations())),
      m_factor_group_index(iter.m_factor_group_index),
      m_translation_index(iter.m_translation_index) {}

PermuteIterator::PermuteIterator(SupercellSymInfo const &_sym_info,
                                 Index _factor_group_index,
                                 Index _translation_index)
    : m_sym_info(&_sym_info),
      m_trans_permute(&(m_sym_info->translation_permutations())),
      m_factor_group_index(_factor_group_index),
      m_translation_index(_translation_index) {}

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
  return m_trans_permute->at(m_translation_index);
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
  if (m_translation_index == m_trans_permute->size()) {
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
    m_translation_index = m_trans_permute->size();
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

void swap(PermuteIterator &a, PermuteIterator &b) {
  std::swap(a.m_sym_info, b.m_sym_info);
  std::swap(a.m_trans_permute, b.m_trans_permute);
  std::swap(a.m_factor_group_index, b.m_factor_group_index);
  std::swap(a.m_translation_index, b.m_translation_index);
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
