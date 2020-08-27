#include "TestOrbits.hh"
#include <algorithm>
#include <iterator>
#include "casm/app/enum.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"

namespace test {

  TestPrimPeriodicIntegralClusterOrbits::TestPrimPeriodicIntegralClusterOrbits(
    const PrimClex &primclex,
    jsonParser _specs) :
    specs(_specs) {

    try {
      // TODO: update with ClusterSpecs
      // Make PrimPeriodicIntegralClusterOrbit
      // make_prim_periodic_orbits(
      //   primclex.shared_prim(),
      //   specs,
      //   alloy_sites_filter,
      //   primclex.crystallography_tol(),
      //   std::back_inserter(orbits),
      //   primclex.log());
    }
    catch(std::exception &e) {
      default_err_log().error("'make_prim_periodic_orbits' failed in TestPrimPeriodicIntegralClusterOrbits");
      //default_err_log() << "prim: \n" << primclex.prim() << std::endl;
      default_err_log() << "specs: \n" << specs << std::endl;
      default_err_log() << e.what() << std::endl;
      throw e;
    }
  }
}
