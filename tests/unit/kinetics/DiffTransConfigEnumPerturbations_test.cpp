#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
//#include "casm/kinetics/DiffTransConfigEnumPerturbations.hh"
//#include "casm/kinetics/DiffTransEnumEquivalents.hh"

/// What is being used to test it:
#include "casm/clex/PrimClex.hh"
#include "casm/app/AppIO_impl.hh"
#include "Common.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/symmetry/InvariantSubgroup_impl.hh"
#include "casm/symmetry/SubOrbits_impl.hh"
#include "casm/kinetics/DiffTransConfigEnumPerturbations.hh"
//#include "casm/casm_io/VaspIO.hh"

using namespace CASM;
using namespace test;

typedef Orbit <
Kinetics::DiffusionTransformation,
         Kinetics::PrimPeriodicDiffTransSymCompare > PrimPeriodicDiffTransOrbit;

// Test 1: PrimPeriodic -> ScelPeriodic orbits
template<typename OrbitType, typename ElementType>
void test_1(
  const Configuration &config,
  const Configuration &prim_config,
  const std::vector<OrbitType> &orbits,
  std::vector<ElementType> &scel_generators,
  std::vector<Index> &orbit_index,
  std::vector<Index> &scel_suborbit_size) {

  // useful SymGroups
  const Structure &prim = config.prim();
  const auto &prim_fg = prim.factor_group();
  const auto &config_scel_fg = config.supercell().factor_group();
  const auto &prim_config_scel_fg = prim_config.supercell().factor_group();

  std::cout << "Symmetry:\n";
  std::cout << "prim_fg.size(): " << prim_fg.size() << std::endl;
  std::cout << "config_scel_fg.size(): " << config_scel_fg.size() << std::endl;
  std::cout << "prim_config_scel_fg.size(): " << prim_config_scel_fg.size() << std::endl;

  std::cout << "\n!!! begin test 1" << std::endl;
  // Split IntegralCluster orbits according to config.supercell().factor_group()
  Index orbit_i = 0;
  for(const auto &orbit : orbits) {
    std::cout << "\n ----------------- \n";
    std::cout << "begin orbit " << orbit_i << "/" << orbits.size() << std::endl;
    std::vector<ElementType> suborbit_generators;
    make_suborbit_generators(
      orbit,
      prim_fg,
      prim_config_scel_fg,
      std::back_inserter(suborbit_generators));

    // Check results:
    std::cout << "Prototype, orbit " << orbit_i << "  prim-orbit.size(): "
              << orbit.size() << std::endl;
    //std::cout << orbit.prototype() << std::endl;
    std::cout << "Sub-orbit generators: " << std::endl;
    Index suborbit_i = 0;
    Index suborbit_size_sum = 0;
    for(const auto &el : suborbit_generators) {
      Orbit<ElementType, PrimPeriodicSymCompare<ElementType>> suborbit(el, prim_config_scel_fg, orbit.sym_compare());
      std::cout << "sub-orbit " << suborbit_i << "/" << suborbit_generators.size()
                << ", size: " << suborbit.size() << ":" << std::endl;
      //std::cout << el << std::endl;
      scel_suborbit_size.push_back(suborbit.size());
      orbit_index.push_back(orbit_i);
      suborbit_size_sum += suborbit.size();
      ++suborbit_i;
    }
    std::cout << "sum: " << suborbit_size_sum
              << "  prim-orbit.size(): "
              << orbit.size() << std::endl;

    BOOST_CHECK_EQUAL(orbit.size(), suborbit_size_sum);
    ++orbit_i;
    std::copy(suborbit_generators.begin(), suborbit_generators.end(), std::back_inserter(scel_generators));
  }

}

// Test 2: ScelPeriodic orbits -> Config orbits
template<typename OrbitType, typename ElementType>
void test_2(
  const Configuration &prim_config,
  const std::vector<OrbitType> &orbits,
  const std::vector<ElementType> &scel_generators,
  std::vector<ElementType> &config_generators,
  std::vector<Index> &orbit_index,
  std::vector<Index> &scel_suborbit_size) {

  const Structure &prim = prim_config.prim();
  std::vector<PermuteIterator> prim_config_fg = prim_config.factor_group();
  SymGroup _prim_config_fg = make_sym_group(prim.lattice(), prim_config_fg);
  ScelPeriodicSymCompare<ElementType> prim_config_scel_sym_compare(
    prim_config.supercell().prim_grid(),
    prim_config.crystallography_tol());

  std::vector<Index> config_suborbit_size;
  std::cout << "\n!!! begin test 2" << std::endl;
  // Split IntegralCluster orbits according to config occupation
  for(Index el_i = 0; el_i < scel_generators.size(); ++el_i) {
    std::cout << "\n ----------------- \n";
    std::cout << "begin scel_generator " << el_i << "/" << scel_generators.size() << std::endl;
    const auto &el = scel_generators[el_i];

    std::vector<ElementType> suborbit_generators;
    make_suborbit_generators(
      el,
      prim_config.supercell(),
      prim_config_fg.begin(),
      prim_config_fg.end(),
      std::back_inserter(suborbit_generators));

    // Check results:
    std::cout << "Generator, " << el_i << "  orbit, " << orbit_index[el_i] << "  scel-orbit.size() / prim-orbit.size(): "
              << scel_suborbit_size[el_i] << "/" << orbits[orbit_index[el_i]].size() << std::endl;
    //std::cout << el << std::endl;
    std::cout << "Sub-orbit generators: " << std::endl;
    Index config_suborbit_i = 0;
    Index config_suborbit_size_sum = 0;
    for(const auto &config_el : suborbit_generators) {
      Orbit<ElementType, ScelPeriodicSymCompare<ElementType>> suborbit(config_el, _prim_config_fg, prim_config_scel_sym_compare);
      std::cout << "config sub-orbit " << config_suborbit_i << "/" << suborbit_generators.size()
                << ", size: " << suborbit.size() << ":" << std::endl;
      //std::cout << config_el << std::endl;
      config_suborbit_size.push_back(suborbit.size());
      config_suborbit_size_sum += suborbit.size();
      ++config_suborbit_i;
    }
    std::cout << "sum: " << config_suborbit_size_sum
              << "  sub-orbit.size() * scel volume: "
              << scel_suborbit_size[el_i]*prim_config.supercell().volume() << std::endl;

    if(el.size() == 0) {
      BOOST_CHECK_EQUAL(1, config_suborbit_size_sum);
    }
    else {
      BOOST_CHECK_EQUAL(scel_suborbit_size[el_i]*prim_config.supercell().volume(), config_suborbit_size_sum);
    }
    std::copy(suborbit_generators.begin(), suborbit_generators.end(), std::back_inserter(config_generators));
  }
  std::cout << "  config_generators.size(): " << config_generators.size() << std::endl;
}

// Test 3: PrimPeriodic -> Config orbits,
//   using make_suborbit_generators_slow and primitive Configuration
template<typename OrbitType, typename ElementType>
void test_3(
  const Configuration &prim_config,
  const std::vector<OrbitType> &orbits,
  const std::vector<ElementType> &config_generators) {

  std::cout << "\n!!! begin test 3" << std::endl;
  std::vector<ElementType> config_generators_slow;
  {
    make_suborbit_generators_slow(
      orbits.begin(),
      orbits.end(),
      prim_config,
      std::back_inserter(config_generators_slow));
    std::cout << "  config_generators_slow.size(): " << config_generators_slow.size() << std::endl;
    BOOST_CHECK_EQUAL(config_generators.size(), config_generators_slow.size());
  }

}

// Test 4: PrimPeriodic -> Config orbits,
//   using make_suborbit_generators and non-primitive Configuration
//   (Need a new test config: current prim_config is lower symmetry than config)
template<typename OrbitType>
void test_4(
  const Configuration &config,
  const std::vector<OrbitType> &orbits) {

  typedef typename OrbitType::Element ElementType;

  std::cout << "\n!!! begin test 4" << std::endl;
  std::vector<ElementType> config_generators_nonprim_slow;
  std::vector<ElementType> config_generators_nonprim;
  {
    make_suborbit_generators_slow(
      orbits.begin(),
      orbits.end(),
      config,
      std::back_inserter(config_generators_nonprim_slow));
    std::cout << "  config_generators_nonprim_slow.size(): " << config_generators_nonprim_slow.size() << std::endl;

    make_suborbit_generators(
      orbits.begin(),
      orbits.end(),
      config,
      std::back_inserter(config_generators_nonprim));
    std::cout << "  config_generators_nonprim.size(): " << config_generators_nonprim.size() << std::endl;

    BOOST_CHECK_EQUAL(config_generators_nonprim_slow.size(), config_generators_nonprim.size());
  }
}



BOOST_AUTO_TEST_SUITE(DiffTransConfigEnumPerturbationsTest)

BOOST_AUTO_TEST_CASE(Test0) {

  test::ZrOProj proj;
  proj.check_init();
  proj.check_composition();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);
  const Structure &prim = primclex.prim();
  const Lattice &lat = prim.lattice();

  fs::path bspecs_path = "tests/unit/kinetics/bspecs_0.json";
  jsonParser bspecs {bspecs_path};

  // Make PrimPeriodicIntegralClusterOrbit
  std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
  make_prim_periodic_orbits(
    primclex.prim(),
    bspecs,
    alloy_sites_filter,
    primclex.crystallography_tol(),
    std::back_inserter(orbits),
    primclex.log());

  print_clust(orbits.begin(), orbits.end(), std::cout, PrototypePrinter<IntegralCluster>());

  // Make PrimPeriodicDiffTransOrbit
  std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
  Kinetics::make_prim_periodic_diff_trans_orbits(
    orbits.begin() + 2,
    orbits.begin() + 4,
    primclex.crystallography_tol(),
    std::back_inserter(diff_trans_orbits));

  print_clust(
    diff_trans_orbits.begin(),
    diff_trans_orbits.end(),
    std::cout,
    PrototypePrinter<Kinetics::DiffusionTransformation>());

  // Make test config
  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();
  Supercell scel {&primclex, Lattice(2 * a, 2 * b, 3 * c)};
  Configuration config(scel);
  config.init_occupation();
  config.init_displacement();
  config.init_deformation();
  config.init_specie_id();
  config.set_occupation({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0});

  // Make test prim_config
  Configuration prim_config = config.primitive().in_canonical_supercell();

  // IntegralCluster tests
  {
    std::vector<IntegralCluster> scel_generators;
    std::vector<Index> orbit_index;
    std::vector<Index> scel_suborbit_size;
    std::vector<IntegralCluster> config_generators;

    Kinetics::DiffusionTransformation diff_trans_prototype = diff_trans_orbits[0].prototype();

    /// Test bubble checkers
    ///Make various supercells
    Supercell scel1 {&primclex, Lattice(2 * a, 2 * b, 3 * c)};
    Supercell scel2 {&primclex, Lattice(4 * a, 4 * b, 3 * c)};
    Supercell scel3 {&primclex, Lattice(8 * a, 2 * b, 3 * c)};
    Supercell scel4 {&primclex, Lattice(4 * a, 3 * b, 4 * c)};
    std::vector<Supercell> scel_list;
    scel_list.push_back(scel1);
    scel_list.push_back(scel2);
    scel_list.push_back(scel3);
    scel_list.push_back(scel4);
    ///make local orbits
    fs::path local_bspecs_path = "tests/unit/kinetics/local_bspecs_0.json";
    jsonParser local_bspecs {local_bspecs_path};
    std::vector<LocalIntegralClusterOrbit> local_orbits;
    make_local_orbits(
      diff_trans_prototype,
      local_bspecs,
      alloy_sites_filter,
      primclex.crystallography_tol(),
      std::back_inserter(local_orbits),
      primclex.log());

    BOOST_CHECK_EQUAL(Kinetics::has_local_bubble_overlap(local_orbits, scel1), 1);
    BOOST_CHECK_EQUAL(Kinetics::has_local_bubble_overlap(local_orbits, scel2), 0);
    BOOST_CHECK_EQUAL(Kinetics::has_local_bubble_overlap(local_orbits, scel3), 1);
    BOOST_CHECK_EQUAL(Kinetics::has_local_bubble_overlap(local_orbits, scel4), 1);
    std::vector<Supercell> result = Kinetics::viable_supercells(local_orbits, scel_list);
    BOOST_CHECK_EQUAL(*(result.begin()) == scel2, 1);

    //test_1(config, prim_config, orbits, scel_generators, orbit_index, scel_suborbit_size);
    //test_2(prim_config, orbits, scel_generators, config_generators, orbit_index, scel_suborbit_size);
    //test_3(prim_config, orbits, config_generators);
    //test_4(config, orbits);

    Configuration config_scel2(scel2);
    config_scel2.init_occupation();
    config_scel2.init_displacement();
    config_scel2.init_deformation();
    config_scel2.init_specie_id();
    config_scel2.set_occupation({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0
                                });


    //FCC TESTING PROCEDURE
    test::FCCTernaryProj fccproj;
    fccproj.check_init();
    fccproj.check_composition();

    Logging fcclogging = Logging::null();
    PrimClex fccprimclex(fccproj.dir, fcclogging);
    const Structure &fccprim = fccprimclex.prim();
    const Lattice &fcclat = fccprim.lattice();

    fs::path fccbspecs_path = "tests/unit/kinetics/bspecs_0.json";
    jsonParser fccbspecs {fccbspecs_path};

    // Make PrimPeriodicIntegralClusterOrbit
    std::vector<PrimPeriodicIntegralClusterOrbit> fccorbits;
    make_prim_periodic_orbits(
      fccprimclex.prim(),
      fccbspecs,
      alloy_sites_filter,
      fccprimclex.crystallography_tol(),
      std::back_inserter(fccorbits),
      fccprimclex.log());

    // Make PrimPeriodicDiffTransOrbit
    std::vector<Kinetics::PrimPeriodicDiffTransOrbit> fccdiff_trans_orbits;
    Kinetics::make_prim_periodic_diff_trans_orbits(
      fccorbits.begin() + 2,
      fccorbits.begin() + 4,
      fccprimclex.crystallography_tol(),
      std::back_inserter(fccdiff_trans_orbits));

    Kinetics::DiffusionTransformation fccdiff_trans_prototype = fccdiff_trans_orbits[4].prototype();
    Eigen::Vector3d a1, b1, c1;
    std::tie(a1, b1, c1) = fccprimclex.prim().lattice().vectors();
    Supercell fccscel {&fccprimclex, Lattice(2 * a1, 2 * b1, 2 * c1)};
    Configuration l12config(fccscel);
    l12config.init_occupation();
    l12config.init_displacement();
    l12config.init_deformation();
    l12config.init_specie_id();
    l12config.set_occupation({0, 0, 0, 1, 1, 0, 0, 0});
    fs::path l12_local_bspecs_path = "tests/unit/kinetics/l12_local_bspecs_0.json";
    jsonParser l12_local_bspecs {l12_local_bspecs_path};

    //In this config there should be 2 options to place the nearest neighbor hop
    // one toward the majority L12 atom and one towards minority L12 atom
    //given a cutoff radius of 5 angstroms and only looking at local point and pair clusters
    //There are the following unique perturbations: (This project still has 3 possible occupants)
    // Hop towards minority L12 surround 5 angst radius with Majority L12
    // Hop towards minority L12 surround 5 angst radius with Minority L12
    // Hop towards minority L12 1/3 of sites around hop with Minority L12 on multiplicity 2 site
    // Hop towards minority L12 1/3 of sites around hop with Minority L12 on one multiplicity 4 site
    // Hop towards minority L12 2/3 of sites around hop with Minority L12 on multiplicity 4 sites
    // Hop towards minority L12 2/3 of sites around hop with Minority L12 on one multiplicity 2 site and one multiplicity 4 site
    //Due to high incidence of periodicity the other orientation of the hop results in the same DiffTransConfigs

    std::set<Kinetics::DiffTransConfiguration> collection;
    Kinetics::DiffTransConfigEnumPerturbations enumerator(l12config, fccdiff_trans_orbits[4], l12_local_bspecs);
    collection.insert(enumerator.begin(), enumerator.end());
    std::cout << collection.size() << std::endl;
    for(auto &dtc : collection) {
      std::cout << "From config" << dtc.sorted().from_config() << std::endl;
      std::cout << "To config " << dtc.sorted().to_config() << std::endl;
    }
  }

  // DiffusionTransformation tests
  {
    std::vector<Kinetics::DiffusionTransformation> scel_generators;
    std::vector<Index> orbit_index;
    std::vector<Index> scel_suborbit_size;
    std::vector<Kinetics::DiffusionTransformation> config_generators;

    //test_1(config, prim_config, diff_trans_orbits, scel_generators, orbit_index, scel_suborbit_size);
    //test_2(prim_config, diff_trans_orbits, scel_generators, config_generators, orbit_index, scel_suborbit_size);
    //test_3(prim_config, diff_trans_orbits, config_generators);
    //test_4(config, diff_trans_orbits);
  }

  BOOST_AUTO_TEST_SUITE_END();
}
