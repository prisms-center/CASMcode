#ifndef CASM_IntegralClusterSymCompareTraits_impl
#define CASM_IntegralClusterSymCompareTraits_impl

#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clusterography/IntegralClusterSymCompareTraits.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/crystallography/SymTools.hh"
// #include "casm/symmetry/Orbit_impl.hh"

namespace CASM {

  template <typename SymCompareType>
  xtal::UnitCellCoord IntegralClusterOrbitTraits<SymCompareType>::position(
    IntegralCluster const &clust,
    SymCompareType const &sym_compare) {
    return clust[0];
  }

  template <typename SymCompareType>
  IntegralCluster IntegralClusterOrbitTraits<SymCompareType>::copy_apply(
    SymOp const &op,
    IntegralCluster const &clust,
    SymCompareType const &sym_compare) {
    return CASM::sym::copy_apply(op, clust, sym_compare.prim());
  }

  template <typename SymCompareType>
  ClusterInvariants IntegralClusterOrbitTraits<SymCompareType>::make_invariants(
    IntegralCluster const &clust,
    SymCompareType const &sym_compare) {
    return ClusterInvariants {clust};
  }

}

#endif
