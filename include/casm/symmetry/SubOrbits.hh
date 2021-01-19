#ifndef CASM_SubOrbits
#define CASM_SubOrbits

namespace CASM {

// class SymGroup;

/// \brief Output the orbit generators necessary to construct the sub-orbits
/// corresponding to group -> subgroup symmetry breaking
template <typename GroupOpIterator, typename SubgroupOpIterator>
class MakeSubOrbitGenerators {
 public:
  MakeSubOrbitGenerators(GroupOpIterator group_begin, GroupOpIterator group_end,
                         SubgroupOpIterator subgroup_begin,
                         SubgroupOpIterator subgroup_end)
      : m_group_begin(group_begin),
        m_group_end(group_end),
        m_subgroup_begin(subgroup_begin),
        m_subgroup_end(subgroup_end) {}

  // /// Output generating elements for the sub-orbits corresponding to group ->
  // subgroup symmetry breaking template<typename Element, typename
  // SymCompareType, typename ElementOutputIterator> ElementOutputIterator
  // operator()(
  //   Element const &element,
  //   SymCompareType const &sym_compare,
  //   ElementOutputIterator result) const;

  /// Output generating elements for the sub-orbits corresponding to group ->
  /// subgroup symmetry breaking
  template <typename Element, typename CopyApplyFunctionType,
            typename PrepareFunctionType, typename InvariantSubgroupOpIterator,
            typename ElementOutputIterator>
  ElementOutputIterator operator()(
      Element const &element, CopyApplyFunctionType copy_apply_f,
      PrepareFunctionType prepare_f,
      InvariantSubgroupOpIterator invariant_subgroup_begin,
      InvariantSubgroupOpIterator invariant_subgroup_end,
      ElementOutputIterator result) const;

 private:
  GroupOpIterator const m_group_begin;
  GroupOpIterator const m_group_end;
  SubgroupOpIterator const m_subgroup_begin;
  SubgroupOpIterator const m_subgroup_end;
};

template <typename GroupOpIterator, typename SubgroupOpIterator>
MakeSubOrbitGenerators<GroupOpIterator, SubgroupOpIterator>
make_suborbit_generators_f(GroupOpIterator group_begin,
                           GroupOpIterator group_end,
                           SubgroupOpIterator subgroup_begin,
                           SubgroupOpIterator subgroup_end) {
  return MakeSubOrbitGenerators<GroupOpIterator, SubgroupOpIterator>{
      group_begin, group_end, subgroup_begin, subgroup_end};
}

template <typename GroupOpIterator, typename SubgroupOpIterator,
          typename ElementIterator, typename SymCompareType,
          typename ElementOutputIterator>
ElementOutputIterator make_suborbit_generators(
    GroupOpIterator group_begin, GroupOpIterator group_end,
    SubgroupOpIterator subgroup_begin, SubgroupOpIterator subgroup_end,
    ElementIterator element_begin, ElementIterator element_end,
    SymCompareType const &sym_compare, ElementOutputIterator result);

}  // namespace CASM

#endif
