#include "casm/container/Permutation.hh"

#include "casm/misc/algorithm.hh"

namespace CASM {

/* PERMUTATION CLASS DEFINITON:

class Permutation {
private:
  /// Permutation array. If m_perm_array[i]==j, then an entry at index 'j'
before permutation
  /// goes to index 'i' after permutation. Put another way, m_perm_array IS the
result of performing
  /// the permutation to the integer list {0,1,2,...,size()-1}
  std::vector<Index> m_perm_array;

public:
  Permutation(const std::vector<Index> &init_perm): m_perm_array(init_perm){};


  Index size() const { return m_perm_array.size();};

  /// Checks that m_perm_array contains values from 0 to m_perm_array.size()-1
and that no value is repeated bool is_perm() const;

  /// Checks whether any indices remain unchanged by permutation
  bool has_fixed_points() const;

  /// Construct permutation that undoes the permutation performed by 'this'
  Permutation inverse() const;

  /// Rearrange 'this' permutation to form an equivalent permutation for
  /// any list that has already been permuted by trans_perm.
  Permutation transformed_by(const Permutation& trans_perm) const;

  /// const access of m_perm_array for doing low-level permutation algebra
  const Index&  operator[](Index i) const { return m_perm_array[i];};

  /// Generate permuted copy of type-T std::vector
  template<typename T>
  std::vector<T> permute(const std::vector<T> &before_array) const;

  /// Generate inversely permuted copy of type-T std::vector
  template<typename T>
  std::vector<T> ipermute(const std::vector<T> &before_array) const;


};
*/

//**************************************************************

/// Checks that m_perm_array contains values from 0 to m_perm_array.size()-1 and
/// that no value is repeated does not depend on definition of permutation
/// convention
bool Permutation::is_perm() const {
  for (Index i = 0; i < size(); i++) {
    if (!valid_index(m_perm_array[i]) || m_perm_array[i] >= size() ||
        find_index(m_perm_array, m_perm_array[i]) < i)
      return false;
  }
  return true;
}

//**************************************************************

Index Permutation::character() const {
  Index result = 0;
  for (Index i = 0; i < size(); ++i) {
    if (m_perm_array[i] == i) ++result;
  }
  return result;
}

//**************************************************************

/// Checks whether any indices remain unchanged by permutation
bool Permutation::has_fixed_points() const {
  for (Index i = 0; i < size(); i++) {
    if (m_perm_array[i] == i) return true;
  }
  return false;
}

//**************************************************************

/// Checks whether any indices remain unchanged by permutation
bool Permutation::is_identity() const {
  for (Index i = 0; i < size(); i++) {
    if (m_perm_array[i] != i) return false;
  }
  return true;
}

//**************************************************************
/// Add new indices that remain unchanged by permutation
void Permutation::append_fixed_points(Index N_new) {
  if (!N_new) return;
  Index orig_size = m_perm_array.size();
  for (Index i = 0; i < N_new; ++i) {
    m_perm_array.push_back(i + orig_size);
  }
}

//**************************************************************

/// Construct permutation that undoes the permutation performed by 'this'
/// Inverse operation is calculated the same, regardless of permutation
/// convention
Permutation Permutation::inverse() const {
  std::vector<Index> im_perm_array(size(), 0);
  for (Index i = 0; i < size(); i++) {
    im_perm_array[m_perm_array[i]] = i;
  }
  return Permutation(std::move(im_perm_array));
}

//**************************************************************
/// Construct permutation of dimension size()*sum(block_dims) that describes the
/// effect of permuting N blocks, where the i'th block has dimension
/// block_dims[i] (N=size()==block_dims.size())
Permutation Permutation::make_block_permutation(
    const std::vector<Index> &block_dims) const {
  assert(block_dims.size() == size());
  std::vector<Index> block_perm;
  block_perm.reserve(sum(block_dims) * size());
  std::vector<Index> i_start(block_dims.size(), 0);
  for (Index i = 0; i + 1 < block_dims.size(); i++) {
    i_start[i + 1] = i_start[i] + block_dims[i];
  }

  for (Index i = 0; i < size(); i++) {
    for (Index j = 0; j < block_dims[m_perm_array[i]]; ++j) {
      block_perm.push_back(i_start[m_perm_array[i]] + j);
    }
  }
  return Permutation(std::move(block_perm));
}

//**************************************************************
/// Given N distinct objects labeled from 0 to N-1,
/// a permutation 'P_permute' that physically permutes the objects (with labels)
/// in terms of their labels, and a permutation 'L_permute' that permutes their
/// labels only, rewrite 'P_permute' in terms of the relabeling induced by
/// 'L_permute' Rearrange 'this' permutation to form an equivalent permutation
/// for any list that has already been permuted by trans_perm. Does not (nearly
/// certain of this) depend on permutation convention
Permutation Permutation::transformed_by(const Permutation &trans_perm) const {
  // Equivalent to Permutation(trans_perm*(*this)*trans_perm.inverse());
  // There's probably a faster element-wise implementation, but it would be
  // confusing
  return Permutation(
      trans_perm.permute(permute(trans_perm.inverse().m_perm_array)));
}

//**************************************************************

//**************************************************************

/// Get the effective Permutation induced by applying RHS, followed by *this
/// I don't think this depends on definition of permutation convention
Permutation Permutation::operator*(const Permutation &RHS) const {
  // not sure if this is obvious.
  return Permutation(permute(RHS.m_perm_array));
}

//**************************************************************

}  // namespace CASM
