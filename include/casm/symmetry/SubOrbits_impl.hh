#ifndef CASM_SubOrbits_impl
#define CASM_SubOrbits_impl

#include "casm/symmetry/InvariantSubgroup_impl.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SubOrbits.hh"

namespace CASM {

  namespace SubOrbits_impl {

    struct OpCompare {
      bool operator()(SymOp const &LHS, SymOp const &RHS) const {
        return LHS.master_group_index() < RHS.master_group_index();
      }
      bool operator()(PermuteIterator const &LHS, PermuteIterator const &RHS) const {
        return LHS < RHS;
      }
    };

    SymOp const &to_sym_op(SymOp const &sym_op) {
      return sym_op;
    }

    SymOp to_sym_op(PermuteIterator const &permute_it) {
      return permute_it.sym_op();
    }

  }

  // --- MakeSubOrbitGenerators ---
  //

  /// Output generating elements for the sub-orbits corresponding to group -> subgroup symmetry breaking
  ///
  /// \param element An orbit generating element w.r.t. [group_begin, group_end)
  /// \param sym_compare SymCompareType functor for comparing elements and applying symmetry
  /// \param result Output iterator outputs sub-orbit generating elements w.r.t. [subgroup_begin, subgroup_end)
  ///
  /// The invariant subgroup of element w.r.t. [group_begin, group_end) is generated as an
  /// temporary intermediate. If it is already available use the other overload.
  template<typename GroupOpIterator, typename SubgroupOpIterator>
  template<typename Element, typename SymCompareType, typename ElementOutputIterator>
  ElementOutputIterator MakeSubOrbitGenerators<GroupOpIterator, SubgroupOpIterator>::operator()(
    Element const &element,
    SymCompareType const &sym_compare,
    ElementOutputIterator result) const {

    typedef typename std::remove_reference<decltype(*std::declval<GroupOpIterator &>())>::type OpType;
    std::vector<OpType> invariant_subgroup;
    make_invariant_subgroup(element, m_group_begin, m_group_end, sym_compare, std::back_inserter(invariant_subgroup));
    return (*this)(element, sym_compare, invariant_subgroup.begin(), invariant_subgroup.end(), result);
  }

  /// Output generating elements for the sub-orbits corresponding to group -> subgroup symmetry breaking
  ///
  /// \param element An orbit generating element w.r.t. [group_begin, group_end)
  /// \param invariant_subgroup_begin, invariant_subgroup_end The subgroup of [group_begin, group_end)
  ///        that leaves element invariant as determined by "sym_compare"
  /// \param sym_compare SymCompareType functor for comparing elements and applying symmetry
  /// \param result Output iterator outputs sub-orbit generating elements w.r.t. [subgroup_begin, subgroup_end)
  ///
  template<typename GroupOpIterator, typename SubgroupOpIterator>
  template<typename Element, typename SymCompareType, typename InvariantSubgroupOpIterator, typename ElementOutputIterator>
  ElementOutputIterator MakeSubOrbitGenerators<GroupOpIterator, SubgroupOpIterator>::operator()(
    Element const &element,
    SymCompareType const &sym_compare,
    InvariantSubgroupOpIterator invariant_subgroup_begin,
    InvariantSubgroupOpIterator invariant_subgroup_end,
    ElementOutputIterator result) const {

    using namespace SubOrbits_impl;

    // Each vector element will contain the set of operations (a subset of [group_begin, group-end)
    // that transform "element" into one of the new sub-orbits:
    //
    //     suborbit_elements = subgroup_op * group_op * invariant_group_op * element
    //
    // The product (group_op * invariant_group_op) generates all the operations that transform
    // element into the equivalent elements. Then, for any particular group_op, the product with
    // subgroup_op generates all the operations that transform element into a particular sub-orbit.
    // ways that element transforms into a particular sub-orbit.
    typedef typename std::remove_reference<decltype(*std::declval<GroupOpIterator &>())>::type OpType;
    std::vector<std::set<OpType, OpCompare>> transform_to_suborbit_sets;

    // For all operations in [group_begin, group_end)
    for(auto group_it = m_group_begin; group_it != m_group_end; ++group_it) {

      // Check if the operation (*group_it) has been added to a set
      bool already_found = false;
      for(auto const &current_set : transform_to_suborbit_sets) {
        if(current_set.count(*group_it)) {
          already_found = true;
        }
      }
      if(already_found) {
        continue;
      }

      // If the operation (*group_it) has not yet been added to a set,
      // use it to find all operations that transform element to the same sub-orbit
      std::set<OpType, OpCompare> next_set;
      auto invariant_subgroup_it = invariant_subgroup_begin;
      for(; invariant_subgroup_it != invariant_subgroup_end; ++invariant_subgroup_it) {
        for(auto subgroup_it = m_subgroup_begin; subgroup_it != m_subgroup_end; ++subgroup_it) {
          next_set.insert((*subgroup_it) * (*group_it) * (*invariant_subgroup_it));
        }
      }
      transform_to_suborbit_sets.push_back(std::move(next_set));
    }

    // Use the first operation in each set to
    // transform element into a sub-orbit generator.
    for(auto const &current_set : transform_to_suborbit_sets) {
      *result++ = sym_compare.prepare(sym_compare.copy_apply(to_sym_op(*current_set.begin()), element));
    }
    return result;
  }

  template<typename GroupOpIterator, typename SubgroupOpIterator, typename ElementIterator, typename SymCompareType, typename ElementOutputIterator>
  ElementOutputIterator make_suborbit_generators(
    GroupOpIterator group_begin,
    GroupOpIterator group_end,
    SubgroupOpIterator subgroup_begin,
    SubgroupOpIterator subgroup_end,
    ElementIterator element_begin,
    ElementIterator element_end,
    SymCompareType const &sym_compare,
    ElementOutputIterator result) {
    MakeSubOrbitGenerators<GroupOpIterator, SubgroupOpIterator> f {
      group_begin, group_end, subgroup_begin, subgroup_end};
    for(auto it = element_begin; it != element_end; ++it) {
      result = f(*it, sym_compare, result);
    }
    return result;
  }

}

#endif
