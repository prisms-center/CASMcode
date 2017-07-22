#ifndef CASM_ScelSubOrbits
#define CASM_ScelSubOrbits

namespace CASM {

  class Supercell;

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to Prim Structure -> Supercell symmetry breaking
  template<typename Element, typename ElementOutputIterator, typename PermuteIteratorIt>
  ElementOutputIterator make_suborbit_generators(
    const Element &element,
    const Supercell &scel,
    PermuteIteratorIt subgroup_begin,
    PermuteIteratorIt subgroup_end,
    ElementOutputIterator result);
}

#endif
