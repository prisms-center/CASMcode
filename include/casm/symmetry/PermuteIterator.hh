#ifndef PermuteIterator_HH
#define PermuteIterator_HH

#include "casm/container/Array.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/container/Permutation.hh"
#include "casm/misc/Comparisons.hh"

namespace CASM {
  class PrimGrid;

  /** \ingroup SymOp
   *  @{
   */

  /// Permutation bidirectional Iterator class
  ///   Can iterate over all combined factor group and translation permutations for a Supercell
  ///
  ///   The sequence:
  ///   for( it = scel.permute_begin(); it != scel.permute_end(); ++it) {
  ///     after_array = it->permute(before_array);
  ///   }
  ///
  ///   replicates this sequence of permutations:
  ///   for( f=0; f<scel->factor_group().size(); f++) {
  ///     for( t=0; t<scel->translation_permute().size(); t++) {
  ///       after_array = scel->translation_permute(t).permute(scel->factor_group_permute(f).permute(before_array));
  ///     }
  ///   }
  ///
  ///   Bidirectional iterators are supposed to be input/output iterators,
  ///     but this is actually only an input iterator... meaning operator* returns by value
  ///
  ///
  class PermuteIterator :
    public std::iterator <std::bidirectional_iterator_tag, PermuteIterator>,
    public Comparisons<PermuteIterator> {

    /// permutation representation of factor group acting on sites of the supercell
    SymGroupRep::RemoteHandle m_fg_permute_rep;

    /// m_prim_grid holds permutation representation of lattice translations acting on sites of the supercell
    PrimGrid const *m_prim_grid;

    /// m_trans_permute points to the Array<Permutation> of translation permutations inside of m_prim_grid (to provide faster access)
    Array<Permutation> const *m_trans_permute;

    Index m_factor_group_index;
    Index m_translation_index;

  public:

    PermuteIterator();

    PermuteIterator(const PermuteIterator &iter);

    PermuteIterator(SymGroupRep::RemoteHandle _fg_permute_rep,
                    const PrimGrid &_prim_grid,
                    Index _factor_group_index,
                    Index _translation_index);

    PermuteIterator &operator=(PermuteIterator iter);

    /// Returns the combination of factor_group permutation and translation permutation
    const PermuteIterator &operator*() const;

    /// Returns the combination of factor_group permutation and translation permutation
    Permutation combined_permute() const;

    /// Apply the combined factor_group permutation and translation permutation being pointed at
    template<typename T>
    ReturnArray<T> permute(const Array<T> &before_array) const;

    /// Return the index into m_factor_group_permute of the factor group op being pointed at
    Index factor_group_index() const;

    /// Return the index into m_prim_grid of the translation being pointed at
    Index translation_index() const;

    /// Return the factor group permutation being pointed at
    const Permutation &factor_group_permute() const;

    /// Return the translation permutation being pointed at
    const Permutation &translation_permute() const;

    /// gets the SymOp for the current operation, defined by translation_op[trans_index]*factor_group_op[fg_index]
    /// i.e, equivalent to application of the factor group operation, FOLLOWED BY application of the translation operation
    SymOp sym_op()const;

    /// Index-wise permutation defined via:
    ///    after_permutation[i] = before_permutation[permutat_iterator.permute_ind(i)];
    Index permute_ind(Index i) const;

    /// Return after_array[i], given i and before_array
    template<typename T>
    const T &permute_by_bit(Index i, const Array<T> &before_array) const;

    bool operator<(const PermuteIterator &iter) const;

    // prefix ++PermuteIterator
    PermuteIterator &operator++();

    // postfix PermuteIterator++
    PermuteIterator operator++(int);

    // prefix --PermuteIterator
    PermuteIterator &operator--();

    // postfix PermuteIterator--
    PermuteIterator operator--(int);

    /// Iterator to next beginning of next factor group operation
    /// skipping all of the intervening operations that differ only by a translation
    PermuteIterator begin_next_fg_op() const;

    PermuteIterator inverse() const;

    PermuteIterator operator*(const PermuteIterator &RHS) const;

    jsonParser &to_json(jsonParser &json) const;
    void from_json(const jsonParser &json);

    friend void swap(PermuteIterator &a, PermuteIterator &b);

  private:

    friend Comparisons<PermuteIterator>;

    bool _eq(const PermuteIterator &iter) const;

  };

  jsonParser &to_json(const PermuteIterator &clust, jsonParser &json);
  void from_json(PermuteIterator &clust, const jsonParser &json);

  /** @} */
}
#endif
