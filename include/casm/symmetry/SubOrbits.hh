#ifndef CASM_SubOrbits
#define CASM_SubOrbits

namespace CASM {

  class SymGroup;

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to group -> subgroup symmetry breaking
  template<typename CopyApplyType>
  class MakeSubOrbitGenerators {
  public:

    MakeSubOrbitGenerators(const SymGroup &group,
                           const SymGroup &subgroup,
                           const CopyApplyType &copy_apply_f):
      m_group(group),
      m_subgroup(subgroup),
      m_copy_apply_f(copy_apply_f) {}

    template<typename OrbitType, typename ElementOutputIterator>
    ElementOutputIterator operator()(
      const OrbitType &orbit,
      ElementOutputIterator result) const;

    template<typename Element, typename SymCompareType, typename ElementOutputIterator>
    ElementOutputIterator operator()(
      const Element &element,
      const SymCompareType &sym_compare,
      ElementOutputIterator result) const;

    template<typename Element, typename ElementOutputIterator>
    ElementOutputIterator operator()(
      const Element &element,
      const SymGroup &invariant_subgroup,
      ElementOutputIterator result) const;

  private:

    const SymGroup &m_group;
    const SymGroup &m_subgroup;
    const CopyApplyType m_copy_apply_f;
  };

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to group -> subgroup symmetry breaking
  template<typename Element, typename CopyApplyElementType, typename ElementOutputIterator>
  ElementOutputIterator make_suborbit_generators(
    const Element &element,
    const SymGroup &invariant_subgroup,
    const SymGroup &group,
    const SymGroup &subgroup,
    const CopyApplyElementType &copy_apply_f,
    ElementOutputIterator result);

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to group -> subgroup symmetry breaking
  template<typename OrbitType, typename CopyApplyElementType, typename ElementOutputIterator>
  ElementOutputIterator make_suborbit_generators(
    const OrbitType &orbit,
    const SymGroup &group,
    const SymGroup &subgroup,
    const CopyApplyElementType &copy_apply_f,
    ElementOutputIterator result);

}

#endif
