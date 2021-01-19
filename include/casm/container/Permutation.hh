// Permutation.hh

#ifndef PERMUTATION_HH
#define PERMUTATION_HH

#include <cassert>
#include <numeric>
#include <ostream>
#include <vector>

#include "casm/global/definitions.hh"

namespace CASM {

/// Permutation is a an operator that permutes indices or items in a list.
/// Permutations are defined such that for an initial array
///        std::vector<THINGS> before;
///
/// A Permutation
///        Permutation my_perm;
///
/// and an output array after permutation
///        std::vector<THINGS> after = my_perm.permute(before);
///
/// The following relation is satisfied
///        after[i] = before[my_perm[i]];
///
/// This is written
///        after = my_perm*before
///
/// Permutations combine
///        after = perm2*perm1*before
///        after = comb_perm*before; comb_perm = perm2*perm1
///
/// Which satisfies
///        after[i] = before[perm1[perm2[i]]]
///        after[i] = before[comb_perm[i]]; comb_perm[i] = perm1[perm2[i]]
///

class Permutation {
 private:
  /// Permutation array. If m_perm_array[i]==j, then an entry at index 'j'
  /// before permutation goes to index 'i' after permutation. Put another way,
  /// m_perm_array IS the result of performing the permutation to the list
  /// {0,1,2,...,size()-1}
  std::vector<Index> m_perm_array;

 public:
  explicit Permutation(Index N) : m_perm_array(N) {
    std::iota(m_perm_array.begin(), m_perm_array.end(), 0);
  }

  template <typename Iterator>
  Permutation(Iterator begin, Iterator end) : m_perm_array(begin, end) {}
  explicit Permutation(const std::vector<Index> &init_perm)
      : m_perm_array(init_perm) {}
  explicit Permutation(std::vector<Index> &&init_perm)
      : m_perm_array(std::move(init_perm)) {}

  Index size() const { return m_perm_array.size(); }

  const std::vector<Index> &perm_array() const { return m_perm_array; }

  /// Checks that m_perm_array contains values from 0 to m_perm_array.size()-1
  /// and that no value is repeated
  bool is_perm() const;

  /// Checks whether any indices remain unchanged by permutation
  bool has_fixed_points() const;

  /// Character of permutation is number of fixed points
  Index character() const;

  /// Checks whether permutation is identity (i.e., m_perm_aray[i]==i for all i)
  bool is_identity() const;

  /// Add new indices that remain unchanged by permutation
  void append_fixed_points(Index N_new);

  /// Construct permutation that undoes the permutation performed by 'this'
  Permutation inverse() const;

  /// Construct permutation of dimension size()*block_dims.sum() that describes
  /// the effect of permuting N blocks, where the i'th block has dimension
  /// block_dims[i] (N=size()==block_dims.size())
  Permutation make_block_permutation(
      const std::vector<Index> &block_dims) const;

  /// Rearrange 'this' permutation to form an equivalent permutation for
  /// any list that has already been permuted by trans_perm.
  Permutation transformed_by(const Permutation &trans_perm) const;

  /// const access of m_perm_array for doing low-level permutation algebra
  const Index &operator[](Index i) const { return m_perm_array[i]; }

  /// non-const access of m_perm_array for doing low-level permutation algebra
  Index &set(Index i) { return m_perm_array[i]; }

  /// Generate permuted copy of indexed container
  template <typename Container>
  Container permute(const Container &_before) const;

  /// Generate inversely permuted copy of type-T std::vector
  template <typename Container>
  Container ipermute(const Container &_before) const;

  template <typename Container>
  Container operator*(const Container &_before) const {
    return this->permute(_before);
  }

  Permutation operator*(const Permutation &RHS) const;
};

std::ostream &operator<<(std::ostream &, const Permutation &);

//************************************************************************************************************************************//
//************************************************************************************************************************************//

/// Generate permuted copy of type-T std::vector
/// THIS IMPLEMENTATION DEPENDS ON DEFINITION OF PERMUTATION CONVENTION
/// Note: To switch conventions, swap names of permute and ipermute
template <typename Container>
Container Permutation::permute(const Container &before_vec) const {
  assert(before_vec.size() == size() &&
         "WARNING: You're trying to permute a Container with an incompatible "
         "permutation!");

  Container after_vec(before_vec);
  // after_vec.reserve(size());

  for (Index i = 0; i < size(); i++) {
    after_vec[i] = before_vec[m_perm_array[i]];
  }
  return after_vec;
}

//**************************************************************

/// Generate inversely permuted copy of type-T std::vector
/// THIS IMPLEMENTATION DEPENDS ON DEFINITION OF PERMUTATION CONVENTION
/// Note: To switch conventions, swap names of permute and ipermute
template <typename Container>
Container Permutation::ipermute(const Container &before_array) const {
  assert(before_array.size() == size() &&
         "WARNING: You're trying to permute an Container with an incompatible "
         "permutation!");

  Container after_array(before_array);
  // after_array.reserve(size());

  for (Index i = 0; i < size(); i++) {
    if (i != m_perm_array[i]) {
      after_array[m_perm_array[i]] = before_array[i];
    }
  }
  return after_array;
}

}  // namespace CASM

#endif
