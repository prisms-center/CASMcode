#include "casm/clusterography/ClusterInvariants.hh"
#include "casm/clusterography/IntegralClusterSymCompareTraits_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clusterography/IntegralCluster.hh"

namespace CASM {

  template class AperiodicSymCompare<IntegralCluster>;
  template class PrimPeriodicSymCompare<IntegralCluster>;
  template class ScelPeriodicSymCompare<IntegralCluster>;
  // template class WithinScelSymCompare<IntegralCluster>;


  // IntegralCluster traits<WithinScelSymCompare<IntegralCluster>>::bring_within(
  //                                                              IntegralCluster clust,
  // WithinScelSymCompare<IntegralCluster> const &sym_compare) {
  //   for(Index i = 0; i < clust.size(); ++i) {
  //     clust[i] = sym_compare.m_bring_within_f(clust[i]);
  //   }
  //   return clust;
  // }
  //
  // IntegralCluster traits<WithinScelSymCompare<IntegralCluster>>::copy_apply(
  //                                                              SymOp const &op,
  //                                                              IntegralCluster const &clust,
  // WithinScelSymCompare<IntegralCluster> const &sym_compare) {
  //   return CASM::sym::copy_apply(op, clust, sym_compare.prim());
  // }
  //
  // WithinScelClusterInvariants traits<WithinScelSymCompare<IntegralCluster>>::make_invariants(
  //                                                                          IntegralCluster clust,
  // WithinScelSymCompare<IntegralCluster> const &sym_compare) {
  //   return WithinScelClusterInvariants {clust, sym_compare.transf_mat()};
  // }
}
