#include "casm/container/Permutation.hh"

namespace CASM {

  /* PERMUTATION CLASS DEFINITON:

  class Permutation {
  private:
    /// Permutation array. If m_perm_array[i]==j, then an entry at index 'j' before permutation
    /// goes to index 'i' after permutation. Put another way, m_perm_array IS the result of performing
    /// the permutation to the integer list {0,1,2,...,size()-1}
    Array<Index> m_perm_array;

  public:
    Permutation(const Array<Index> &init_perm): m_perm_array(init_perm){};


    Index size() const { return m_perm_array.size();};

    /// Checks that m_perm_array contains values from 0 to m_perm_array.size()-1 and that no value is repeated
    bool is_perm() const;

    /// Checks whether any indices remain unchanged by permutation
    bool has_fixed_points() const;

    /// Construct permutation that undoes the permutation performed by 'this'
    Permutation inverse() const;

    /// Rearrange 'this' permutation to form an equivalent permutation for
    /// any list that has already been permuted by trans_perm.
    Permutation transformed_by(const Permutation& trans_perm) const;

    /// const access of m_perm_array for doing low-level permutation algebra
    const Index&  operator[](Index i) const { return m_perm_array[i];};

    /// Generate permuted copy of type-T Array
    template<typename T>
    Array<T> permute(const Array<T> &before_array) const;

    /// Generate inversely permuted copy of type-T Array
    template<typename T>
    Array<T> ipermute(const Array<T> &before_array) const;


  };
  */

  //**************************************************************

  /// Checks that m_perm_array contains values from 0 to m_perm_array.size()-1 and that no value is repeated
  /// does not depend on definition of permutation convention
  bool Permutation::is_perm() const {
    for(Index i = 0; i < size(); i++) {
      if(!valid_index(m_perm_array[i]) || m_perm_array[i] >= size() || m_perm_array.find(m_perm_array[i]) < i)
        return false;
    }
    return true;
  }


  //**************************************************************

  /// Checks whether any indices remain unchanged by permutation
  bool Permutation::has_fixed_points() const {
    for(Index i = 0; i < size(); i++) {
      if(m_perm_array[i] == i)
        return true;
    }
    return false;
  }

  //**************************************************************
  /// Add new indices that remain unchanged by permutation
  void Permutation::append_fixed_points(Index N_new) {
    if(!N_new)
      return;
    m_perm_array.append(Array<Index>::sequence(m_perm_array.size(), m_perm_array.size() + N_new - 1));
  }

  //**************************************************************

  /// Construct permutation that undoes the permutation performed by 'this'
  /// Inverse operation is calculated the same, regardless of permutation convention
  Permutation Permutation::inverse() const {
    Array<Index> im_perm_array(size(), 0);
    for(Index i = 0; i < size(); i++) {
      im_perm_array[m_perm_array[i]] = i;
    }
    return Permutation(ReturnArray<Index>(im_perm_array));
  }

  //**************************************************************
  /// Construct permutation of dimension size()*block_dims.sum() that describes the effect of permuting
  /// N blocks, where the i'th block has dimension block_dims[i] (N=size()==block_dims.size())
  Permutation Permutation::make_block_permutation(const Array<Index> &block_dims)const {
    assert(block_dims.size() == size());
    Array<Index> block_perm;
    block_perm.reserve(block_dims.sum()*size());
    Array<Index> i_start(block_dims.size(), 0);
    for(Index i = 0; i + 1 < block_dims.size(); i++) {
      i_start[i + 1] = i_start[i] + block_dims[i];
    }

    for(Index i = 0; i < size(); i++) {
      block_perm.append(Array<Index>::sequence(i_start[m_perm_array[i]], i_start[m_perm_array[i]] + block_dims[m_perm_array[i]] - 1));
    }
    return Permutation(ReturnArray<Index>(block_perm));
  }

  //**************************************************************
  /// Given N distinct objects labeled from 0 to N-1,
  /// a permutation 'P_permute' that physically permutes the objects (with labels) in terms of their labels,
  /// and a permutation 'L_permute' that permutes their labels only,
  /// rewrite 'P_permute' in terms of the relabeling induced by 'L_permute'
  /// Rearrange 'this' permutation to form an equivalent permutation for
  /// any list that has already been permuted by trans_perm.
  /// Does not (nearly certain of this) depend on permutation convention
  Permutation Permutation::transformed_by(const Permutation &trans_perm) const {

    // Equivalent to Permutation(trans_perm*(*this)*trans_perm.inverse());
    // There's probably a faster element-wise implementation, but it would be confusing
    return Permutation(ReturnArray<Index>(trans_perm.permute(permute(trans_perm.inverse().m_perm_array))));
  }

  //**************************************************************

  //**************************************************************

  /// Get the effective Permutation induced by applying RHS, followed by *this
  /// I don't think this depends on definition of permutation convention
  Permutation Permutation::operator*(const Permutation &RHS) const {
    //not sure if this is obvious.
    return Permutation(ReturnArray<Index>(permute(RHS.m_perm_array)));
  }

  //**************************************************************

  jsonParser &Permutation::to_json(jsonParser &json) const {
    return CASM::to_json(m_perm_array, json);
    /*    json.put_array();
    for(Index i = 0; i < size(); i++)
      json.push_back(m_perm_array[i]);
    */
    //return json;
  }

  //**************************************************************

  void Permutation::from_json(const jsonParser &json) {
    try {

      CASM::from_json(m_perm_array, json);
      /*m_perm_array.resize(json.size());
      for(Index i = 0; i < json.size(); i++)
      from_json(m_perm_array[i], json[i]);*/
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

  //**************************************************************

  std::ostream &operator<<(std::ostream &out, const Permutation &perm) {
    out << perm.perm_array();
    return out;
  }

  //**************************************************************

  jsonParser &to_json(const Permutation &perm, jsonParser &json) {
    return perm.to_json(json);
  }

  //**************************************************************

  void from_json(Permutation &perm, const jsonParser &json) {
    try {
      perm.from_json(json);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }
}

