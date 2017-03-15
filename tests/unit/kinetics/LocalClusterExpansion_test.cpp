#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"

/// What is being used to test it:
#include "casm/clex/PrimClex.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/AppIO_impl.hh"
#include "Common.hh"
#include "casm/completer/Handlers.hh"
#include "casm/app/casm_functions.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffusionTransformationEnum.hh"
#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"

using namespace CASM;
using namespace test;

BOOST_AUTO_TEST_SUITE(LocalClusterExpansionTest)

BOOST_AUTO_TEST_CASE(Test0) {

  test::ZrOProj proj;
  proj.check_init();
  proj.check_composition();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);

  fs::path bspecs_path = "tests/unit/kinetics/bspecs_0.json";
  jsonParser bspecs {bspecs_path};

  std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
  make_prim_periodic_orbits(
    primclex.prim(),
    bspecs,
    alloy_sites_filter,
    primclex.crystallography_tol(),
    std::back_inserter(orbits),
    primclex.log());

  std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
  Kinetics::make_prim_periodic_diff_trans_orbits(orbits.begin() + 4, orbits.begin() + 7, primclex.crystallography_tol(), std::back_inserter(diff_trans_orbits));
  Kinetics::DiffusionTransformation trans = diff_trans_orbits[0].prototype();

  fs::path local_bspecs_path = "tests/unit/kinetics/local_bspecs_0.json";
  jsonParser local_bspecs {local_bspecs_path};

  std::vector<LocalIntegralClusterOrbit> local_orbits;
  make_local_orbits(
    trans,
    local_bspecs,
    alloy_sites_filter,
    primclex.crystallography_tol(),
    std::back_inserter(local_orbits),
    primclex.log());


  std::cout << trans << std::endl;
  PrototypePrinter<IntegralCluster> printer;
  print_clust(local_orbits.begin(), local_orbits.end(), std::cout, printer);



}

BOOST_AUTO_TEST_SUITE_END()
