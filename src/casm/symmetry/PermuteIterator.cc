#include "casm/symmetry/PermuteIterator.hh"

#include "casm/crystallography/PrimGrid.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/clex/Supercell.hh"
#include "casm/casm_io/jsonParser.hh"

namespace CASM {

  PermuteIterator::PermuteIterator() {}

  PermuteIterator::PermuteIterator(const PermuteIterator &iter) :
    m_sym_info(iter.m_sym_info),
    m_trans_permute(&(m_sym_info->prim_grid().translation_permutations())),
    m_factor_group_index(iter.m_factor_group_index),
    m_translation_index(iter.m_translation_index) {

  }

  PermuteIterator::PermuteIterator(SupercellSymInfo const &_sym_info,
                                   Index _factor_group_index,
                                   Index _translation_index) :
    m_sym_info(&_sym_info),
    m_trans_permute(&(m_sym_info->prim_grid().translation_permutations())),
    m_factor_group_index(_factor_group_index),
    m_translation_index(_translation_index) {
  }

  const PrimGrid &PermuteIterator::prim_grid() const {
    return sym_info().prim_grid();
  }

  PermuteIterator &PermuteIterator::operator=(PermuteIterator iter) {
    swap(*this, iter);
    return *this;
  }

  /// Returns a reference to this
  const PermuteIterator &PermuteIterator::operator*() const {
    return *this;
  }

  /// Returns a pointer to this
  const PermuteIterator *PermuteIterator::operator->() const {
    return this;
  }

  /// Returns the combination of factor_group permutation and translation permutation
  Permutation PermuteIterator::combined_permute() const {
    return translation_permute() * factor_group_permute();
  }

  /// Return the index into Supercell::factor_group_permute of the factor group op being pointed at
  SupercellSymInfo const &PermuteIterator::sym_info() const {
    return *m_sym_info;
  }

  /// Return the index into Supercell::factor_group_permute of the factor group op being pointed at
  SymGroup const &PermuteIterator::factor_group() const {
    return sym_info().factor_group();
  }

  /// Return the index into Supercell::factor_group_permute of the factor group op being pointed at
  Index PermuteIterator::factor_group_index() const {
    return m_factor_group_index;
  }

  /// Return the index into Supercell::translation_permute of the translation being pointed at
  Index PermuteIterator::translation_index() const {
    return m_translation_index;
  }

  /// Return the factor group permutation being pointed at
  const Permutation &PermuteIterator::factor_group_permute() const {
    return *(sym_info().site_permutation_symrep()[m_factor_group_index]->permutation());
  }

  /// Return the translation permutation being pointed at
  const Permutation &PermuteIterator::translation_permute() const {
    return m_trans_permute->at(m_translation_index);
  }

  /// Returns representation of current operation corresponding to species transformation on sublattice b
  SymOpRepresentation const &PermuteIterator::occ_rep(Index b) const {
    return *sym_info().occ_symreps()[b][factor_group_index()];
  }

  /// Returns representation of current operation
  /// corresponding to local DoF specified by _key on sublattice b
  SymOpRepresentation const &PermuteIterator::local_dof_rep(DoFKey const &_key, Index b) const {
    return *sym_info().local_dof_symreps(_key)[b][factor_group_index()];
  }

  /// Returns representation of current operation corresponding to global DoF specified by _key
  SymOpRepresentation const &PermuteIterator::global_dof_rep(DoFKey const &_key) const {
    return *sym_info().global_dof_symrep(_key)[factor_group_index()];
  }

  SymOp PermuteIterator::sym_op()const {
    return prim_grid().sym_op(m_translation_index) * sym_info().factor_group()[m_factor_group_index];
  }

  Index PermuteIterator::permute_ind(Index i) const {
    return factor_group_permute()[translation_permute()[i]];
  }


  bool PermuteIterator::operator<(const PermuteIterator &iter) const {
    if(this->factor_group_index() == iter.factor_group_index()) {
      return this->translation_index() < iter.translation_index();
    }
    return this->factor_group_index() < iter.factor_group_index();
  }

  bool PermuteIterator::eq_impl(const PermuteIterator &iter) const {
    if(m_sym_info == iter.m_sym_info &&
       m_factor_group_index == iter.m_factor_group_index &&
       m_translation_index == iter.m_translation_index)
      return true;
    return false;
  }

  // prefix ++PermuteIterator
  PermuteIterator &PermuteIterator::operator++() {
    m_translation_index++;
    if(m_translation_index == m_trans_permute->size()) {
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
    if(m_translation_index == 0) {
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

  // Might be able to do this in a faster way, but would be much harder to understand
  PermuteIterator PermuteIterator::inverse() const {
    PermuteIterator it(*this);
    // Finding the inverse factor_group operation is straightforward
    it.m_factor_group_index = factor_group().ind_inverse(factor_group_index());

    // Easiest way to get the new translation is just to compare the tau of the
    // inverse of the 'total' sym_op (described by *this), to the inverse of the
    // untranslated symop (described by m_fg_site_permutation_symrep.sym_op(it.m_factor_group_index))
    // Result is the portion of the inverse sym_op that needs to be described by a prim_grid translation
    it.m_translation_index =
      prim_grid().find_cart(sym_op().inverse().tau() - factor_group()[it.factor_group_index()].tau());

    return it;
  }

  PermuteIterator PermuteIterator::operator*(const PermuteIterator &RHS) const {
    PermuteIterator it(*this);
    // Finding the inverse factor_group operation is straightforward
    it.m_factor_group_index = factor_group().ind_prod(factor_group_index(), RHS.factor_group_index());

    // Easiest way to get the new translation is just to compare the tau of the
    // 'total' sym_op (described by (*this).sym_op()*RHS.sym_op()), to the
    // untranslated symop product (described by m_fg_site_permutation_symrep.sym_op(it.factor_group_index()))
    // Result is the portion of the product sym_op that needs to be described by a prim_grid translation
    it.m_translation_index =
      prim_grid().find_cart((sym_op() * RHS.sym_op()).tau() - factor_group()[it.factor_group_index()].tau());

    return it;
  }

  std::ostream &operator<<(std::ostream &sout, const PermuteIterator &op) {
    sout << "(" << op.factor_group_index() << ", " << op.prim_grid().unitcell(op.translation_index()).transpose() << ")";
    return sout;
  }

  /// \brief Returns a SymGroup generated from a range of PermuteIterator
  ///
  /// \param begin,end A range of PermuteIterator
  ///
  /// - The result is sorted
  /// - The result uses the Supercell lattice for periodic comparisons
  template<typename PermuteIteratorIt>
  SymGroup make_point_group(PermuteIteratorIt begin, PermuteIteratorIt end) {
    SymGroup result;
    result.set_lattice(begin->prim_grid().scel_lattice());
    for(; begin != end; ++begin) {
      Index f = begin->sym_op().index();
      Index i;
      for(i = 0; i < result.size(); ++i) {
        if(f == result[i].index())
          break;
      }
      if(i == result.size()) {
        result.push_back((begin->sym_op()).no_trans());
        std::cout << "pushed back op " << result.back().index() << "\n";
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
  template<typename PermuteIteratorIt>
  SymGroup make_sym_group(PermuteIteratorIt begin, PermuteIteratorIt end) {
    SymGroup sym_group;
    sym_group.set_lattice(begin->prim_grid().scel_lattice());
    while(begin != end) {
      sym_group.push_back(begin->sym_op());
      ++begin;
    }
    sym_group.sort();
    return sym_group;
  }

  template SymGroup make_sym_group(
    PermuteIterator begin,
    PermuteIterator end);
  template SymGroup make_sym_group(
    std::vector<PermuteIterator>::const_iterator begin,
    std::vector<PermuteIterator>::const_iterator end);
  template SymGroup make_sym_group(
    std::vector<PermuteIterator>::iterator begin,
    std::vector<PermuteIterator>::iterator end);

  template SymGroup make_point_group(
    PermuteIterator begin,
    PermuteIterator end);
  template SymGroup make_point_group(
    std::vector<PermuteIterator>::const_iterator begin,
    std::vector<PermuteIterator>::const_iterator end);
  template SymGroup make_point_group(
    std::vector<PermuteIterator>::iterator begin,
    std::vector<PermuteIterator>::iterator end);

  void swap(PermuteIterator &a, PermuteIterator &b) {
    std::swap(a.m_sym_info, b.m_sym_info);
    std::swap(a.m_trans_permute, b.m_trans_permute);
    std::swap(a.m_factor_group_index, b.m_factor_group_index);
    std::swap(a.m_translation_index, b.m_translation_index);
  }

  jsonParser &to_json(const PermuteIterator &it, jsonParser &json) {
    json.put_obj();
    json["factgrp"] = it.factor_group_index();
    json["trans"] = it.translation_index();
    return json;
  }

  PermuteIterator jsonConstructor<PermuteIterator>::from_json(
    const jsonParser &json,
    const SupercellSymInfo &scel_info) {
    return scel_info.permute_it(json["factgrp"].get<Index>(), json["trans"].get<Index>());
  }

}
