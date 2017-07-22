#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/symmetry/SubOrbits.hh"
#include "casm/symmetry/SubOrbits_impl.hh"

#include "casm/symmetry/ScelSubOrbits.hh"
#include "casm/symmetry/ScelSubOrbits_impl.hh"

#include "casm/symmetry/ConfigSubOrbits.hh"
#include "casm/symmetry/ConfigSubOrbits_impl.hh"

/// What is being used to test it:
#include "Common.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"


using namespace CASM;

BOOST_AUTO_TEST_SUITE(SubOrbitsTest)

/// Test individual parts of MakeConfigSubOrbitGenerators
BOOST_AUTO_TEST_CASE(Test0) {

  /// Make test project
  BOOST_CHECK_EQUAL(true, true);
  test::ZrOProj proj;
  proj.check_init();
  proj.check_composition();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);
  const Structure &prim = primclex.prim();
  const Lattice &lat = prim.lattice();
  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();
  BOOST_CHECK_EQUAL(true, true);


  // Make PrimPeriodicIntegralClusterOrbit
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
  BOOST_CHECK_EQUAL(true, true);
  BOOST_CHECK_EQUAL(orbits.size(), 74);

  //print_clust(orbits.begin(), orbits.end(), std::cout, PrototypePrinter<IntegralCluster>());

  // Make PrimPeriodicDiffTransOrbit
  std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
  Kinetics::make_prim_periodic_diff_trans_orbits(
    orbits.begin() + 2,
    orbits.begin() + 4,
    primclex.crystallography_tol(),
    std::back_inserter(diff_trans_orbits));
  BOOST_CHECK_EQUAL(true, true);
  BOOST_CHECK_EQUAL(diff_trans_orbits.size(), 4);

  /*
  std::cout << "diff_trans_orbits: \n" << std::endl;
  print_clust(
    diff_trans_orbits.begin(),
    diff_trans_orbits.end(),
    std::cout,
    PrototypePrinter<Kinetics::DiffusionTransformation>());
  */

  Printer<Kinetics::DiffusionTransformation> diff_trans_printer;
  typedef std::vector<Kinetics::DiffusionTransformation> DiffTransVec;
  typedef ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> ScelDiffTransSymCompare;

  // Test 0 Step-by-step of internals of MakeConfigSubOrbitGenerators
  {
    /*
    Supercell scel {&primclex, Lattice(2 * a, 1 * b, 1 * c)};
    Configuration config(scel);
    config.set_occupation({0, 0, 0, 0, 1, 1, 0, 0});
    */

    Supercell scel {&primclex, Lattice(3 * a, 2 * b, 1 * c)};
    Configuration config(scel);
    config.set_occupation({
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      1, 1, 1, 1, 1, 1,
      0, 0, 0, 0, 0, 0
    });
    ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> scel_sym_compare(
      config.supercell().prim_grid(),
      config.crystallography_tol());


    Configuration prim_config = config.primitive();
    std::vector<PermuteIterator> prim_config_fg = prim_config.factor_group();
    ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> prim_sym_compare(
      prim_config.supercell().prim_grid(),
      prim_config.crystallography_tol());

    BOOST_CHECK_EQUAL(true, true);

    //std::cout << "config: \n" << config << std::endl;
    //std::cout << "prim_config: \n" << prim_config << std::endl;

    auto print_res = [&](std::string name, const DiffTransVec & diff_trans_vec, ScelDiffTransSymCompare sym_compare) {
      std::cout << name << ":" << std::endl;
      Index index = 0;
      for(const auto &el : diff_trans_vec) {
        std::cout << "index: " << index << std::endl;
        diff_trans_printer.print(sym_compare.prepare(el), std::cout);
        ++index;
      }
      std::cout << std::endl;
    };

    /// prim -> prim_config.supercell() symmetry breaking
    DiffTransVec prim_scel_suborbit_generators;
    {
      MakeSubOrbitGenerators gen(
        prim_config.prim().factor_group(),
        prim_config.supercell().factor_group());
      for(const auto &orbit : diff_trans_orbits) {
        gen(orbit, std::back_inserter(prim_scel_suborbit_generators));
      }
      BOOST_CHECK_EQUAL(true, true);

      //print_res("prim_scel_suborbit_generators", prim_scel_suborbit_generators, prim_sym_compare);
    }

    /// prim_config.supercell() -> prim_config symmetry breaking
    std::vector<Kinetics::DiffusionTransformation> prim_config_suborbit_generators;
    {
      for(const auto &el : prim_scel_suborbit_generators) {
        make_suborbit_generators(
          el,
          prim_config.supercell(),
          prim_config_fg.begin(),
          prim_config_fg.end(),
          std::back_inserter(prim_config_suborbit_generators));
      }
      BOOST_CHECK_EQUAL(true, true);

      //print_res("prim_config_suborbit_generators", prim_config_suborbit_generators, prim_sym_compare);
    }

    /// prim_config.supercell() -> config.supercell() symmetry breaking
    std::vector<Kinetics::DiffusionTransformation> config_suborbit_generators;
    {

      //std::cout << "prim_config.supercell().factor_group().size(): "
      //  << prim_config.supercell().factor_group().size() << std::endl;
      //std::cout << "config.supercell().factor_group().size(): "
      //  << config.supercell().factor_group().size() << std::endl;

      MakeSubOrbitGenerators gen(
        prim_config.supercell().factor_group(),
        config.supercell().factor_group());
      for(const auto &el : prim_config_suborbit_generators) {
        gen(el, prim_sym_compare, std::back_inserter(config_suborbit_generators));
      }
      BOOST_CHECK_EQUAL(true, true);
      BOOST_CHECK_EQUAL(config_suborbit_generators.size(), 14);

      //print_res("config_suborbit_generators 0", config_suborbit_generators, scel_sym_compare);
    }

    // test all combined in MakeConfigSubOrbitGenerators
    {
      std::vector<Kinetics::DiffusionTransformation> config_suborbit_generators;
      MakeConfigSubOrbitGenerators{config}(
        diff_trans_orbits.begin(),
        diff_trans_orbits.end(),
        std::back_inserter(config_suborbit_generators));
      BOOST_CHECK_EQUAL(true, true);
      BOOST_CHECK_EQUAL(config_suborbit_generators.size(), 14);

      //print_res("config_suborbit_generators 1", config_suborbit_generators, scel_sym_compare);
    }

    // test all combined in MakeConfigSubOrbitGenerators
    {
      std::vector<Kinetics::DiffusionTransformation> config_suborbit_generators;
      make_suborbit_generators(
        diff_trans_orbits.begin(),
        diff_trans_orbits.end(),
        config,
        std::back_inserter(config_suborbit_generators));
      BOOST_CHECK_EQUAL(true, true);
      BOOST_CHECK_EQUAL(config_suborbit_generators.size(), 14);

      //print_res("config_suborbit_generators 2", config_suborbit_generators, scel_sym_compare);
    }

  }
}

BOOST_AUTO_TEST_SUITE_END()
