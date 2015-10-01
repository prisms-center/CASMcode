// Permutation.hh

#ifndef PERMUTATION_HH
#define PERMUTATION_HH

#include "casm/container/Array.hh"

namespace CASM {

  /// Permutation is a an operator that permutes indices or items in a list.
  /// Permutations are defined such that for an initial array
  ///        Array<THINGS> before;
  ///
  /// A Permutation
  ///        Permutation my_perm;
  ///
  /// and an output array after permutation
  ///        Array<THINGS> after = my_perm.permute(before);
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
    /// Permutation array. If m_perm_array[i]==j, then an entry at index 'j' before permutation
    /// goes to index 'i' after permutation. Put another way, m_perm_array IS the result of performing
    /// the permutation to the list {0,1,2,...,size()-1}
    Array<Index> m_perm_array;

  public:
    Permutation(Index N): m_perm_array(Array<Index>::sequence(0, N - 1)) {};

    Permutation(const std::vector<Index> &init_perm) {
      m_perm_array.reserve(init_perm.size());
      for(auto it = init_perm.cbegin(); it != init_perm.cend(); ++it)
        m_perm_array.push_back(*it);
    };
    Permutation(const Array<Index> &init_perm): m_perm_array(init_perm) {};
    Permutation(ReturnArray<Index> &init_perm): m_perm_array(init_perm) {};

    Index size() const {
      return m_perm_array.size();
    };

    const Array<Index> &perm_array() const {
      return m_perm_array;
    };

    /// Checks that m_perm_array contains values from 0 to m_perm_array.size()-1 and that no value is repeated
    bool is_perm() const;

    /// Checks whether any indices remain unchanged by permutation
    bool has_fixed_points() const;

    /// Add new indices that remain unchanged by permutation
    void append_fixed_points(Index N_new);

    /// Construct permutation that undoes the permutation performed by 'this'
    Permutation inverse() const;

    /// Construct permutation of dimension size()*block_dims.sum() that describes the effect of permuting
    /// N blocks, where the i'th block has dimension block_dims[i] (N=size()==block_dims.size())
    Permutation make_block_permutation(const Array<Index> &block_dims)const;

    /// Rearrange 'this' permutation to form an equivalent permutation for
    /// any list that has already been permuted by trans_perm.
    Permutation transformed_by(const Permutation &trans_perm) const;

    /// const access of m_perm_array for doing low-level permutation algebra
    const Index  &operator[](Index i) const {
      return m_perm_array[i];
    };

    /// const access of m_perm_array for doing low-level permutation algebra
    const Index  &at(Index i) const {
      return m_perm_array[i];
    };

    /// Generate permuted copy of type-T Array
    template<typename T>
    ReturnArray<T> permute(const Array<T> &before_array) const;

    /// Generate permuted copy of type-T vector
    template<typename T>
    std::vector<T> permute(const std::vector<T> &before_array) const;

    /// Generate inversely permuted copy of type-T Array
    template<typename T>
    ReturnArray<T> ipermute(const Array<T> &before_array) const;

    template<typename T>
    ReturnArray<T> operator*(const Array<T> &before_array) const {
      return this->permute(before_array);
    }

    Permutation operator*(const Permutation &RHS) const;

    jsonParser &to_json(jsonParser &json) const;
    void from_json(const jsonParser &json);
  };

  jsonParser &to_json(const Permutation &value, jsonParser &json);
  void from_json(Permutation &value, const jsonParser &json);

  std::ostream &operator<<(std::ostream &, const Permutation &);

  //************************************************************************************************************************************//
  //************************************************************************************************************************************//

  /// Generate permuted copy of type-T Array
  /// THIS IMPLEMENTATION DEPENDS ON DEFINITION OF PERMUTATION CONVENTION
  /// Note: To switch conventions, swap names of permute and ipermute
  template<typename T>
  ReturnArray<T> Permutation::permute(const Array<T> &before_array) const {
    assert(before_array.size() == size() && "WARNING: You're trying to permute an Array with an incompatible permutation!");

    Array<T> after_array;
    after_array.reserve(size());

    for(Index i = 0; i < size(); i++) {
      after_array.push_back(before_array[m_perm_array[i]]);
    }
    return after_array;
  }


  //************************************************************************************************************************************//

  /// Generate permuted copy of type-T std::vector
  /// THIS IMPLEMENTATION DEPENDS ON DEFINITION OF PERMUTATION CONVENTION
  /// Note: To switch conventions, swap names of permute and ipermute
  template<typename T>
  std::vector<T> Permutation::permute(const std::vector<T> &before_vec) const {
    assert(before_vec.size() == size() && "WARNING: You're trying to permute an Array with an incompatible permutation!");

    std::vector<T> after_vec;
    after_vec.reserve(size());

    for(Index i = 0; i < size(); i++) {
      after_vec.push_back(before_vec[m_perm_array[i]]);
    }
    return after_vec;
  }

  //**************************************************************

  /// Generate inversely permuted copy of type-T Array
  /// THIS IMPLEMENTATION DEPENDS ON DEFINITION OF PERMUTATION CONVENTION
  /// Note: To switch conventions, swap names of permute and ipermute
  template<typename T>
  ReturnArray<T> Permutation::ipermute(const Array<T> &before_array) const {
    assert(before_array.size() == size() && "WARNING: You're trying to permute an Array with an incompatible permutation!");

    Array<T> after_array(before_array);
    //after_array.reserve(size());

    for(Index i = 0; i < size(); i++) {
      if(i != m_perm_array[i]) {
        after_array[m_perm_array[i]] = before_array[i];
      }
    }
    return after_array;

  }

}


#endif
