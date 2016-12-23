#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/kinetics/DiffusionTransformation.hh"

/// What is being used to test it:
#include "casm/clex/PrimClex.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/app/AppIO.hh"
#include "Common.hh"

using namespace CASM;
using namespace test;

BOOST_AUTO_TEST_SUITE(DiffusionTransformationTest)

BOOST_AUTO_TEST_CASE(Test0) {

  test::ZrOProj proj;
  make_project(proj);
  proj.check_init();
  proj.check_composition();

  //Logging logging = Logging::null();
  Logging logging;
  std::cout << "here 1" << std::endl;
  PrimClex primclex(proj.dir, logging);

  std::cout << "here 1" << std::endl;
  fs::path bspecs_path = "tests/unit/kinetics/bspecs_0.json";
  jsonParser bspecs {bspecs_path};
  std::cout << "here 2: \n" << bspecs << std::endl;

  std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
  make_prim_periodic_orbits(
    primclex.prim(),
    bspecs,
    alloy_sites_filter,
    primclex.crystallography_tol(),
    std::back_inserter(orbits),
    primclex.log());

  std::cout << "#orbits: " << orbits.size() << std::endl;

  print_clust(orbits.begin(), orbits.end(), primclex.log(), ProtoSitesPrinter());

  rm_project(proj);

}

BOOST_AUTO_TEST_SUITE_END()
