#include "casm/symmetry/PermuteIterator.hh"

#include "casm/crystallography/PrimGrid.hh"

namespace CASM {

  PermuteIterator::PermuteIterator() {}

  PermuteIterator::PermuteIterator(const PermuteIterator &iter) :
    m_fg_permute_rep(iter.m_fg_permute_rep),
    m_prim_grid(iter.m_prim_grid),
    m_trans_permute(&(m_prim_grid->translation_permutations())),
    m_factor_group_index(iter.m_factor_group_index),
    m_translation_index(iter.m_translation_index) {

  }

  PermuteIterator::PermuteIterator(SymGroupRep::RemoteHandle _fg_permute_rep,
                                   const PrimGrid &_prim_grid,
                                   Index _factor_group_index,
                                   Index _translation_index) :
    m_fg_permute_rep(_fg_permute_rep),
    m_prim_grid(&_prim_grid),
    m_trans_permute(&(m_prim_grid->translation_permutations())),
    m_factor_group_index(_factor_group_index),
    m_translation_index(_translation_index) {
  }

  PermuteIterator &PermuteIterator::operator=(PermuteIterator iter) {
    swap(*this, iter);
    return *this;
  }

  /// Returns a copy
  const PermuteIterator &PermuteIterator::operator*() const {
    return *this;
  }

  /// Returns the combination of factor_group permutation and translation permutation
  Permutation PermuteIterator::combined_permute() const {
    return translation_permute() * factor_group_permute();
  }

  /// Apply the combined factor_group permutation and translation permutation being pointed at
  template<typename T>
  ReturnArray<T> PermuteIterator::permute(const Array<T> &before_array) const {
    assert(before_array.size() == factor_group_permute().size() && "WARNING: You're trying to permute an Array with an incompatible permutation!");

    Array<T> after_array;
    after_array.reserve(before_array.size());

    for(int i = 0; i < before_array.size(); i++) {
      after_array.push_back(permute_by_bit(i, before_array));
    }
    return after_array;
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
    return *(m_fg_permute_rep[m_factor_group_index]->get_permutation());
  }

  /// Return the translation permutation being pointed at
  const Permutation &PermuteIterator::translation_permute() const {
    return m_trans_permute->at(m_translation_index);
  }

  SymOp PermuteIterator::sym_op()const {
    return (*m_prim_grid).sym_op(m_translation_index) * m_fg_permute_rep.sym_op(m_factor_group_index);
  }

  Index PermuteIterator::permute_ind(Index i) const {
    return factor_group_permute()[ translation_permute()[i] ];
  }

  /// Return after_array[i], given i and before_array
  template<typename T>
  const T &PermuteIterator::permute_by_bit(Index i, const Array<T> &before_array) const {
    return before_array[ factor_group_permute()[ translation_permute()[i] ] ];
  }

  bool PermuteIterator::operator<(const PermuteIterator &iter) const {
    if(this->factor_group_index() == iter.factor_group_index()) {
      return this->translation_index() < iter.translation_index();
    }
    return this->factor_group_index() < iter.factor_group_index();
  }

  bool PermuteIterator::_eq(const PermuteIterator &iter) const {
    if(m_fg_permute_rep == iter.m_fg_permute_rep &&
       m_prim_grid == iter.m_prim_grid &&
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
    it.m_factor_group_index = m_fg_permute_rep.ind_inverse(factor_group_index());

    // Easiest way to get the new translation is just to compare the tau of the
    // inverse of the 'total' sym_op (described by *this), to the inverse of the
    // untranslated symop (described by m_fg_permute_rep.sym_op(it.m_factor_group_index))
    // Result is the portion of the inverse sym_op that needs to be described by a prim_grid translation
    it.m_translation_index =
      m_prim_grid->find_cart(sym_op().inverse().tau() - m_fg_permute_rep.sym_op(it.factor_group_index()).tau());

    return it;
  }

  PermuteIterator PermuteIterator::operator*(const PermuteIterator &RHS) const {
    PermuteIterator it(*this);
    // Finding the inverse factor_group operation is straightforward
    it.m_factor_group_index = m_fg_permute_rep.ind_prod(factor_group_index(), RHS.factor_group_index());

    // Easiest way to get the new translation is just to compare the tau of the
    // 'total' sym_op (described by (*this).sym_op()*RHS.sym_op()), to the
    // untranslated symop product (described by m_fg_permute_rep.sym_op(it.factor_group_index()))
    // Result is the portion of the product sym_op that needs to be described by a prim_grid translation
    it.m_translation_index =
      m_prim_grid->find_cart((sym_op() * RHS.sym_op()).tau() - m_fg_permute_rep.sym_op(it.factor_group_index()).tau());

    return it;
  }

  jsonParser &PermuteIterator::to_json(jsonParser &json) const {
    json.put_obj();
    json["factgrp"] = m_factor_group_index;
    json["trans"] = m_translation_index;
    return json;
  }

  void PermuteIterator::from_json(const jsonParser &json) {
    CASM::from_json(m_factor_group_index, json["factgrp"]);
    CASM::from_json(m_translation_index, json["trans"]);
  }

  void swap(PermuteIterator &a, PermuteIterator &b) {
    std::swap(a.m_fg_permute_rep, b.m_fg_permute_rep);
    std::swap(a.m_prim_grid, b.m_prim_grid);
    std::swap(a.m_trans_permute, b.m_trans_permute);
    std::swap(a.m_factor_group_index, b.m_factor_group_index);
    std::swap(a.m_translation_index, b.m_translation_index);
  }

  jsonParser &to_json(const PermuteIterator &it, jsonParser &json) {
    return it.to_json(json);
  }

  void from_json(PermuteIterator &it, const jsonParser &json) {
    it.from_json(json);
  }

}
