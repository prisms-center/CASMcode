#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/kinetics/DiffTransConfiguration.hh"

/// What is being used to test it:
#include "casm/clex/PrimClex.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/AppIO_impl.hh"
#include "Common.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffusionTransformationEnum.hh"
#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"
#include "casm/clusterography/ClusterOrbits.hh"

using namespace CASM;
using namespace test;

BOOST_AUTO_TEST_SUITE(DiffTransConfigurationTest)

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

  //print_clust(orbits.begin(), orbits.end(), primclex.log(), ProtoSitesPrinter());
  std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
  Kinetics::make_prim_periodic_diff_trans_orbits(orbits.begin()+4, orbits.begin()+5, primclex.crystallography_tol(), std::back_inserter(diff_trans_orbits));
  Kinetics::DiffusionTransformation trans = diff_trans_orbits[0].prototype();
  
  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  Supercell scel {&primclex, Lattice(2*a, 2*b, 2*c)};

  Configuration config(scel);
  config.init_occupation();
  config.init_displacement();
  config.init_deformation();
  config.set_occupation({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1});
  //test Constructor/field accessors
  Kinetics::DiffTransConfiguration dtc(config,trans);
  BOOST_CHECK_EQUAL(dtc.from_config(), config);
  BOOST_CHECK_EQUAL(dtc.diff_trans(), trans);
  std::cout << config << std::endl;
  Configuration tmp {config};

  tmp = dtc.diff_trans().apply_to(tmp);
  std::cout << "WHERE" << std::endl;
  BOOST_CHECK_EQUAL(dtc.to_config(), tmp);
  std::cout << "WHERE" << std::endl;
  //check sorting
  BOOST_CHECK_EQUAL(dtc.is_sorted(), dtc.from_config() < dtc.to_config());
  BOOST_CHECK_EQUAL(dtc.is_sorted(),dtc == dtc.sorted());


  }

BOOST_AUTO_TEST_SUITE_END()
