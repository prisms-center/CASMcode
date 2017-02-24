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
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/kinetics/DoFTransformation.hh"
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
  Kinetics::make_prim_periodic_diff_trans_orbits(orbits.begin() + 4, orbits.begin() + 7, primclex.crystallography_tol(), std::back_inserter(diff_trans_orbits));
  Kinetics::DiffusionTransformation trans = diff_trans_orbits[0].prototype();
  Kinetics::DiffusionTransformation trans2 = diff_trans_orbits[2].prototype();

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  Supercell scel {&primclex, Lattice(2 * a, 2 * b, 2 * c)};

  Configuration config(scel);
  config.init_occupation();
  config.init_displacement();
  config.init_deformation();
  config.init_specie_id();
  //hardcoded occupation for trans to occur is there a way to do this generally?
  config.set_occupation({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1});

  Configuration config2(scel);
  config2.init_occupation();
  config2.init_displacement();
  config2.init_deformation();
  config2.init_specie_id();
  //hardcoded occupation for trans to occur is there a way to do this generally?
  config2.set_occupation({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0});




  //test Constructor/field accessors
  Kinetics::DiffTransConfiguration dtc(config, trans);
  BOOST_CHECK_EQUAL(dtc.from_config(), config);
  BOOST_CHECK_EQUAL(dtc.diff_trans(), trans);
  Configuration tmp {config};

  tmp = dtc.diff_trans().apply_to(tmp);
  BOOST_CHECK_EQUAL(dtc.to_config(), tmp);

  //check comparison
  Kinetics::DiffTransConfiguration dtc2(config2, trans2);
  Kinetics::DiffTransConfiguration dtc3(config2, trans);

  //config > config2 but trans < trans2
  //comparing transformation takes priority
  BOOST_CHECK_EQUAL(dtc < dtc2, trans < trans2);
  //If the above test fails the preference of priority may have changed
  BOOST_CHECK_EQUAL(dtc < dtc3, config < config2);

  //check apply sym
  PermuteIterator it = config.supercell().permute_begin();
  BOOST_CHECK_EQUAL(copy_apply(it, dtc) == dtc, 1);

  ++it;
  ++it;
  ++it;
  ++it;
  ++it;
  ++it;
  ++it;
  ++it;
  Configuration new_config = copy_apply(it, config);
  Kinetics::ScelPeriodicDiffTransSymCompare symcompare(config.supercell().prim_grid(),
                                                       config.supercell().crystallography_tol());

  Kinetics::DiffusionTransformation new_trans =
    symcompare.prepare(copy_apply(it.sym_op(), trans));

  Kinetics::DiffTransConfiguration newdtc(new_config, new_trans);

  BOOST_CHECK_EQUAL(copy_apply(it, dtc) == newdtc, 1);

  //check sorting
  BOOST_CHECK_EQUAL(dtc.is_sorted(), dtc.from_config() < dtc.to_config());
  BOOST_CHECK_EQUAL(dtc.is_sorted(), dtc == dtc.sorted());

  //check canonical form
  BOOST_CHECK_EQUAL(dtc.is_canonical(), 0);
  std::cout << (dtc < dtc.canonical_form()) << std::endl;
  BOOST_CHECK_EQUAL(!dtc.is_canonical(), dtc < dtc.canonical_form());
  BOOST_CHECK_EQUAL(dtc.is_canonical(), dtc == dtc.canonical_form());
  BOOST_CHECK_EQUAL(1, dtc.canonical_form().is_canonical());
  BOOST_CHECK_EQUAL(1, copy_apply(dtc.to_canonical(), dtc) == dtc.canonical_form());
  BOOST_CHECK_EQUAL(1, copy_apply(dtc.from_canonical(), dtc.canonical_form()) == dtc);


}

BOOST_AUTO_TEST_SUITE_END()
