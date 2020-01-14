#ifndef PermuteIterator_HH
#define PermuteIterator_HH

#include <vector>
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SupercellSymInfo.hh"
#include "casm/container/Permutation.hh"
#include "casm/misc/Comparisons.hh"

namespace CASM {

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
    public Comparisons<CRTPBase<PermuteIterator>> {

    SupercellSymInfo const *m_sym_info;

    /// m_trans_permute points to the vector<Permutation> of translation permutations inside of m_sym_info (to provide faster access)
    std::vector<Permutation> const *m_trans_permute;

    Index m_factor_group_index;
    Index m_translation_index;

  public:

    //TODO: What does this even mean? You're asking for a segmentation fault
    PermuteIterator();

    PermuteIterator(const PermuteIterator &iter);

    PermuteIterator(SupercellSymInfo const &_sym_info,
                    Index _factor_group_index,
                    Index _translation_index);


    PermuteIterator &operator=(PermuteIterator iter);

    /// Returns a reference to this -- allows PermuteIterator to be treated as an iterator to PermuteIterator object
    const PermuteIterator &operator*() const;

    /// Returns a pointer to this -- allows PermuteIterator to be treated as an iterator to PermuteIterator object
    const PermuteIterator *operator->() const;

    /// Returns the combination of factor_group permutation and translation permutation
    Permutation combined_permute() const;

    SymGroup const &factor_group() const;

    //TODO: Get rid of this?
    SupercellSymInfo const &sym_info() const;

    /// Return the index into m_factor_group_permute of the factor group op being pointed at
    Index factor_group_index() const;

    /// Return the index into m_prim_grid of the translation being pointed at
    Index translation_index() const;

    /// Return the factor group permutation being pointed at
    const Permutation &factor_group_permute() const;

    /// Return the translation permutation being pointed at
    const Permutation &translation_permute() const;

    /// Returns representation of current operation corresponding to species transformation on sublattice b
    SymOpRepresentation const &occ_rep(Index b) const;

    /// Returns representation of current operation
    /// corresponding to local DoF specified by _key on sublattice b
    SymOpRepresentation const &local_dof_rep(DoFKey const &_key, Index b) const;

    /// Returns representation of current operation corresponding to global DoF specified by _key
    SymOpRepresentation const &global_dof_rep(DoFKey const &_key) const;

    /// gets the SymOp for the current operation, defined by translation_op[trans_index]*factor_group_op[fg_index]
    /// i.e, equivalent to application of the factor group operation, FOLLOWED BY application of the translation operation
    SymOp sym_op()const;

    /// Index-wise permutation defined via:
    ///    after_permutation[i] = before_permutation[permutat_iterator.permute_ind(i)];
    Index permute_ind(Index i) const;

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

    friend void swap(PermuteIterator &a, PermuteIterator &b);

  private:

    friend Comparisons<CRTPBase<PermuteIterator>>;

    bool eq_impl(const PermuteIterator &iter) const;

  };

  /// \brief Output PermuteIterator as (fg_index, i, j, k)
  /* std::ostream &operator<<(std::ostream &sout, const PermuteIterator &op); */

  /// Iterator to next beginning of next factor group operation
  /// skipping all of the intervening operations that differ only by a translation
  template<typename IterType>
  IterType begin_next_fg_op(IterType it, IterType end) {
    Index tfg = it->factor_group_index();
    for(++it; it != end; ++it) {
      if(it->factor_group_index() != tfg)
        break;
    }
    return it;
  }

  /// \brief Returns a SymGroup generated from a container of PermuteIterator
  ///
  /// \param container A container of PermuteIterator
  ///
  /// - The result is sorted
  template<typename PermuteIteratorContainer>
  SymGroup make_point_group(const PermuteIteratorContainer &container, const Lattice &supercell_lattice) {
    return make_point_group(container.begin(), container.end(), supercell_lattice);
  }

  /// \brief Returns a SymGroup generated from a range of PermuteIterator
  template<typename PermuteIteratorIt>
  SymGroup make_point_group(PermuteIteratorIt begin, PermuteIteratorIt end, const Lattice &supercell_lattice);

  /// \brief Returns a SymGroup generated from a container of PermuteIterator
  ///
  /// \param container A container of PermuteIterator
  ///
  /// - The result is sorted
  template<typename PermuteIteratorContainer>
  SymGroup make_sym_group(const PermuteIteratorContainer &container, const Lattice &supercell_lattice) {
    return make_sym_group(container.begin(), container.end(), supercell_lattice);
  }

  /// \brief Returns a SymGroup generated from a range of PermuteIterator
  template<typename PermuteIteratorIt>
  SymGroup make_sym_group(PermuteIteratorIt begin, PermuteIteratorIt end, const Lattice &supercell_lattice);

  jsonParser &to_json(const PermuteIterator &clust, jsonParser &json);


  class SupercellSymInfo;

  template<typename T> struct jsonConstructor;

  template<>
  struct jsonConstructor<PermuteIterator> {
    static PermuteIterator from_json(const jsonParser &json, const SupercellSymInfo &scel_info);
  };

  /** @} */
}
#endif
