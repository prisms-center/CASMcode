#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/symmetry/PermuteIterator.hh"

namespace CASM {

  IntegralClusterOrbitGenerator::IntegralClusterOrbitGenerator(
    IntegralCluster const &_prototype,
    bool _include_subclusters):
    prototype(_prototype), include_subclusters(_include_subclusters) {}

}
