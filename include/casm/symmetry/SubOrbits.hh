#ifndef CASM_SubOrbits
#define CASM_SubOrbits

namespace CASM {

  class SymGroup;
  class Supercell;
  class Configuration;

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

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to group -> subgroup symmetry breaking
  template<typename Element, typename ElementOutputIterator, typename PermuteIteratorIt>
  ElementOutputIterator make_suborbit_generators(
    const Element &element,
    const Supercell &scel,
    PermuteIteratorIt subgroup_begin,
    PermuteIteratorIt subgroup_end,
    ElementOutputIterator result);

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to Prim Structure -> Configuration symmetry breaking
  template<typename OrbitIterator, typename ElementOutputIterator>
  ElementOutputIterator make_suborbit_generators_slow(
    OrbitIterator begin,
    OrbitIterator end,
    const Configuration &config,
    ElementOutputIterator result);

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to Prim Structure -> Configuration symmetry breaking
  template<typename OrbitIterator, typename ElementOutputIterator>
  ElementOutputIterator make_suborbit_generators(
    OrbitIterator begin,
    OrbitIterator end,
    const Configuration &config,
    ElementOutputIterator result);

}

#endif
