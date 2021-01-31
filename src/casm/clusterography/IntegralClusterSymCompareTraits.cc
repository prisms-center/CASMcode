#include "casm/clusterography/ClusterInvariants.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clusterography/IntegralClusterSymCompareTraits_impl.hh"

namespace CASM {

template class AperiodicSymCompare<IntegralCluster>;
template class PrimPeriodicSymCompare<IntegralCluster>;
template class ScelPeriodicSymCompare<IntegralCluster>;

}  // namespace CASM
