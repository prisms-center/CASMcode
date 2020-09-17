#ifndef CASM_ClusterSymCompareDecl
#define CASM_ClusterSymCompareDecl

#include "casm/symmetry/OrbitDecl.hh"

namespace CASM {

  enum class CLUSTER_PERIODICITY_TYPE : int {
    PRIM_PERIODIC = 0, APERIODIC = 1, LOCAL = 1, WITHIN_SCEL = 2, SCEL_PERIODIC = 3
  };

  template<typename T> struct traits;

  /// \brief Template class to be specialized for comparisons with aperiodic symmetry
  template<typename Element, typename U = void>
  class AperiodicSymCompare;

  template<typename Element, typename U = void>
  using LocalSymCompare = AperiodicSymCompare<Element, U>;

  /// \brief Template class to be specialized for comparisons with periodic symmetry
  /// of the primitive lattice
  template<typename Element, typename U = void>
  class PrimPeriodicSymCompare;

  /// \brief Template class to be specialized for comparisons with periodic symmetry
  /// of the supercell lattice
  template<typename Element, typename U = void>
  class ScelPeriodicSymCompare;

  // /// \brief Template class to be specialized for comparisons with periodic symmetry
  // /// of the supercell lattice, for clusters that do not extend outside the supercell
  // template<typename Element, typename U = void>
  // class WithinScelSymCompare;

  template<typename Element>
  using AperiodicOrbit = Orbit<AperiodicSymCompare<Element>>;

  template<typename Element>
  using LocalOrbit = AperiodicOrbit<Element>;

  template<typename Element>
  using ScelPeriodicOrbit = Orbit<ScelPeriodicSymCompare<Element>>;

  template<typename Element>
  using PrimPeriodicOrbit = Orbit<PrimPeriodicSymCompare<Element>>;

  // template<typename Element>
  // using WithinScelOrbit = Orbit<WithinScelSymCompare<Element>>;

}

#endif
