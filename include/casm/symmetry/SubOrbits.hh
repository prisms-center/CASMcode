#ifndef CASM_SubOrbits
#define CASM_SubOrbits

namespace CASM {

  class SymGroup;

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to group -> subgroup symmetry breaking
  class MakeSubOrbitGenerators {
  public:

    MakeSubOrbitGenerators(const SymGroup &group,
                           const SymGroup &subgroup);

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
  };

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to group -> subgroup symmetry breaking
  template<typename Element, typename ElementOutputIterator>
  ElementOutputIterator make_suborbit_generators(
    const Element &element,
    const SymGroup &invariant_subgroup,
    const SymGroup &group,
    const SymGroup &subgroup,
    ElementOutputIterator result);

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to group -> subgroup symmetry breaking
  template<typename Element, typename SymCompareType, typename ElementOutputIterator>
  ElementOutputIterator make_suborbit_generators(
    const Element &element,
    const SymCompareType &sym_compare,
    const SymGroup &group,
    const SymGroup &subgroup,
    ElementOutputIterator result);

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to group -> subgroup symmetry breaking
  template<typename OrbitType, typename ElementOutputIterator>
  ElementOutputIterator make_suborbit_generators(
    const OrbitType &orbit,
    const SymGroup &group,
    const SymGroup &subgroup,
    ElementOutputIterator result);

}

#endif
