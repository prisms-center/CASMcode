#ifndef PermuteIterator_HH
#define PermuteIterator_HH

#include <vector>

#include "casm/container/Permutation.hh"
#include "casm/misc/Comparisons.hh"
#include "casm/symmetry/SupercellSymInfo.hh"
#include "casm/symmetry/SymGroupRep.hh"

namespace CASM {

/** \ingroup SymOp
 *  @{
 */

/// Iterate over all combined factor group and translation permutations for a
/// Supercell
///
/// The sequence:
/// \code
/// SupercellSymInfo sym_info = ...
/// for(it = sym_info.permute_begin(); it != sym_info.permute_end(); ++it) {
///   after_array = it->permute(before_array);
/// }
/// \endcode
///
/// Replicates this sequence of permutations:
/// \code
/// SupercellSymInfo sym_info = ...
/// for( f=0; f<sym_info.factor_group().size(); f++) {
///   for( t=0; t<sym_info.translation_permutations().size(); t++) {
///     after_array =
///     sym_info.translation_permute(t).permute(sym_info.factor_group_permute(f).permute(before_array));
///   }
/// }
/// \endcode
///
/// Bidirectional iterators are supposed to be input/output iterators, but this
/// is actually only an input iterator (meaning operator* returns by value).
///
///
class PermuteIterator
    : public std::iterator<std::bidirectional_iterator_tag, PermuteIterator>,
      public Comparisons<CRTPBase<PermuteIterator>> {
  SupercellSymInfo const *m_sym_info;

  /// Pointer to the vector<Permutation> of translation permutations
  /// inside of m_sym_info (to provide faster access)
  std::vector<Permutation> const *m_trans_permute;

  Index m_factor_group_index;
  Index m_translation_index;

 public:
  PermuteIterator();

  PermuteIterator(const PermuteIterator &iter);

  PermuteIterator(SupercellSymInfo const &_sym_info, Index _factor_group_index,
                  Index _translation_index);

  PermuteIterator &operator=(PermuteIterator iter);

  /// Returns a reference to this -- allows PermuteIterator to be treated as an
  /// iterator to PermuteIterator object
  const PermuteIterator &operator*() const;

  /// Returns a pointer to this -- allows PermuteIterator to be treated as an
  /// iterator to PermuteIterator object
  const PermuteIterator *operator->() const;

  /// Returns the combination of factor_group permutation and translation
  /// permutation
  Permutation combined_permute() const;

  /// Reference the `SupercellSymInfo::factor_group()`
  SymGroup const &factor_group() const;

  /// Reference the SupercellSymInfo containing the operations being pointed at
  SupercellSymInfo const &sym_info() const;

  /// Return the supercell factor group index
  Index factor_group_index() const;

  /// Return the prim factor group index
  Index prim_factor_group_index() const;

  /// Return the index into the supercell translation permutations
  Index translation_index() const;

  /// Return the factor group permutation being pointed at
  const Permutation &factor_group_permute() const;

  /// Return the translation permutation being pointed at
  const Permutation &translation_permute() const;

  /// Returns representation of current operation corresponding to species
  /// transformation on sublattice b
  SymOpRepresentation const &occ_rep(Index b) const;

  /// Returns representation of current operation
  /// corresponding to local DoF specified by _key on sublattice b
  SymOpRepresentation const &local_dof_rep(DoFKey const &_key, Index b) const;

  /// Check if local DoF representation is empty
  bool local_dof_rep_empty(DoFKey const &_key, Index b) const;

  /// Returns representation of current operation corresponding to global DoF
  /// specified by _key
  SymOpRepresentation const &global_dof_rep(DoFKey const &_key) const;

  /// gets the SymOp for the current operation, defined by
  /// translation_op[trans_index]*factor_group_op[fg_index] i.e, equivalent to
  /// application of the factor group operation, FOLLOWED BY application of the
  /// translation operation
  SymOp sym_op() const;

  /// Index-wise permutation defined via:
  ///    after_permutation[i] =
  ///    before_permutation[permute_iterator.permute_ind(i)];
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
  /// skipping all of the intervening operations that differ only by a
  /// translation
  PermuteIterator begin_next_fg_op() const;

  PermuteIterator inverse() const;

  /// Return the PermuteIterator product (equivalent to first applying RHS and
  /// then applying *this).
  PermuteIterator operator*(const PermuteIterator &RHS) const;

  friend void swap(PermuteIterator &a, PermuteIterator &b);

  /// Returns true if the two PermuteIterators share the same instance of
  /// sym_info. This is meant to be used if you want to check that a particular
  /// range of PermuteIterators are referring to the same supercell, to ensure
  /// your iterating over a consistent set of permutations.
  bool is_compatible(const PermuteIterator &other_permute_iterator) const {
    return (&this->sym_info() == &other_permute_iterator->sym_info());
  }

 private:
  friend Comparisons<CRTPBase<PermuteIterator>>;

  bool eq_impl(const PermuteIterator &iter) const;
};

/// Iterator to next beginning of next factor group operation
/// skipping all of the intervening operations that differ only by a translation
template <typename IterType>
IterType begin_next_fg_op(IterType it, IterType end) {
  Index tfg = it->factor_group_index();
  for (++it; it != end; ++it) {
    if (it->factor_group_index() != tfg) break;
  }
  return it;
}

/// \brief Returns a SymGroup generated from a container of PermuteIterator
///
/// \param container A container of PermuteIterator
///
/// - The result is sorted
template <typename PermuteIteratorContainer>
SymGroup make_point_group(const PermuteIteratorContainer &container,
                          const Lattice &supercell_lattice) {
  return make_point_group(container.begin(), container.end(),
                          supercell_lattice);
}

/// \brief Returns a SymGroup generated from a range of PermuteIterator
template <typename PermuteIteratorIt>
SymGroup make_point_group(PermuteIteratorIt begin, PermuteIteratorIt end,
                          const Lattice &supercell_lattice);

/// \brief Returns a SymGroup generated from a container of PermuteIterator
///
/// \param container A container of PermuteIterator
///
/// - The result is sorted
template <typename PermuteIteratorContainer>
SymGroup make_sym_group(const PermuteIteratorContainer &container,
                        const Lattice &supercell_lattice) {
  return make_sym_group(container.begin(), container.end(), supercell_lattice);
}

/// \brief Returns a SymGroup generated from a range of PermuteIterator
template <typename PermuteIteratorIt>
SymGroup make_sym_group(PermuteIteratorIt begin, PermuteIteratorIt end,
                        const Lattice &supercell_lattice);

/// \brief Returns a std::unique_ptr<SymGroup> generated from a container of
/// PermuteIterator
///
/// \param container A container of PermuteIterator
///
/// - The result is sorted
template <typename PermuteIteratorContainer>
std::unique_ptr<SymGroup> make_unique_sym_group(
    const PermuteIteratorContainer &container,
    const Lattice &supercell_lattice) {
  return make_unique_sym_group(container.begin(), container.end(),
                               supercell_lattice);
}

/// \brief Returns a SymGroup generated from a range of PermuteIterator
template <typename PermuteIteratorIt>
std::unique_ptr<SymGroup> make_unique_sym_group(
    PermuteIteratorIt begin, PermuteIteratorIt end,
    const Lattice &supercell_lattice);

/// \brief Make permute group for local property symmetry in a supercell
std::vector<PermuteIterator> make_local_permute_group(
    SymGroup const &local_generating_group, SymGroup const &factor_group,
    SupercellSymInfo const &supercell_sym_info);

/// \brief Filter PermuteIterator to keep only operations consistent
///     with a factor group in a sub-supercell
template <typename PermuteIteratorIt>
std::vector<PermuteIterator> make_allowed_permute(
    PermuteIteratorIt supercell_permute_begin,
    PermuteIteratorIt supercell_permute_end,
    std::set<PermuteIterator> const &subsupercell_factor_group);

/// Return true if the permutation does not mix given sites and other sites
bool site_indices_are_invariant(PermuteIterator const &permute_it,
                                std::set<Index> const &site_indices);

jsonParser &to_json(const PermuteIterator &clust, jsonParser &json);

class SupercellSymInfo;

template <typename T>
struct jsonConstructor;

template <>
struct jsonConstructor<PermuteIterator> {
  static PermuteIterator from_json(const jsonParser &json,
                                   const SupercellSymInfo &scel_info);
};

namespace adapter {

template <typename ToType, typename FromType>
struct Adapter;

/// Convert CASM::PermuteIterator -> CASM::SymOp
template <>
struct Adapter<SymOp, PermuteIterator> {
  SymOp operator()(PermuteIterator const &adaptable) const;
};
}  // namespace adapter

/** @} */
}  // namespace CASM
#endif
