#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/kinetics/DiffTransConfigEnumOccPerturbations.hh"

/// What is being used to test it:
#include "casm/clex/PrimClex.hh"
#include "Common.hh"
#include "casm/app/AppIO_impl.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/kinetics/DiffusionTransformation_impl.hh"
#include "casm/kinetics/DiffTransConfiguration_impl.hh"
#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"

using namespace CASM;
using namespace test;

typedef Orbit <
Kinetics::DiffusionTransformation,
         Kinetics::PrimPeriodicDiffTransSymCompare > PrimPeriodicDiffTransOrbit;

BOOST_AUTO_TEST_SUITE(DiffTransConfigEnumOccPerturbationsTest)

BOOST_AUTO_TEST_CASE(NeighborhoodOverlapTest) {

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
    std::back_inserter(diff_trans_orbits),
    &primclex);
  BOOST_CHECK_EQUAL(true, true);
  BOOST_CHECK_EQUAL(diff_trans_orbits.size(), 4);

  Kinetics::DiffusionTransformation diff_trans_prototype = diff_trans_orbits[0].prototype();

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

  BOOST_CHECK_EQUAL(Kinetics::has_local_bubble_overlap(local_orbits, scel1), 1);
  BOOST_CHECK_EQUAL(Kinetics::has_local_bubble_overlap(local_orbits, scel2), 0);
  BOOST_CHECK_EQUAL(Kinetics::has_local_bubble_overlap(local_orbits, scel3), 1);
  BOOST_CHECK_EQUAL(Kinetics::has_local_bubble_overlap(local_orbits, scel4), 1);
  std::vector<Supercell> result = Kinetics::viable_supercells(local_orbits, scel_list);
  BOOST_CHECK_EQUAL(*(result.begin()) == scel2, 1);

}

BOOST_AUTO_TEST_CASE(ZrOTest) {

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
    std::back_inserter(diff_trans_orbits),
    &primclex);
  BOOST_CHECK_EQUAL(true, true);
  BOOST_CHECK_EQUAL(diff_trans_orbits.size(), 4);

  /*
  print_clust(
    diff_trans_orbits.begin(),
    diff_trans_orbits.end(),
    std::cout,
    PrototypePrinter<Kinetics::DiffusionTransformation>());
  */

  // Make background config
  Supercell _scel {&primclex, Lattice(1 * a, 1 * b, 1 * c)};
  Configuration _config(_scel);
  _config.set_occupation({0, 0, 1, 0});

  //std::cout << "construct background_config" << std::endl;
  Supercell background_scel {&primclex, Lattice(3 * a, 3 * b, 3 * c)};
  Configuration background_config = _config.
                                    fill_supercell(background_scel, primclex.prim().factor_group()).
                                    in_canonical_supercell();
  BOOST_CHECK_EQUAL(true, true);

  /// Construct enumerator
  //std::cout << "construct enumerator" << std::endl;
  fs::path local_bspecs_path = "tests/unit/kinetics/ZrO_local_bspecs_0.json";
  jsonParser local_bspecs {local_bspecs_path};
  Kinetics::DiffTransConfigEnumOccPerturbations enumerator(
    background_config,
    diff_trans_orbits[0],
    local_bspecs);
  BOOST_CHECK_EQUAL(true, true);
  BOOST_CHECK_EQUAL(enumerator.valid(), true);
  BOOST_CHECK_EQUAL(enumerator.step(), 0);
  BOOST_CHECK_EQUAL((enumerator.begin() != enumerator.end()), true);

  /// Enumerate perturbations (may be duplicates at this point)
  //std::cout << "enumerate" << std::endl;
  std::vector<Kinetics::DiffTransConfiguration> collection;
  Index index = 0;
  for(auto it = enumerator.begin(); it != enumerator.end(); ++it) {
    //std::cout << "OUTPUT: " << index << std::endl;
    //std::cout << "  diff_trans: \n" << it->diff_trans() << std::endl;
    //std::cout << "  occ: " << it->from_config().occupation() << std::endl;
    //std::cout << "  occ: " << it->to_config().occupation() << std::endl;
    collection.push_back(*it);
    BOOST_CHECK_EQUAL(it->has_valid_from_occ(), true);
    BOOST_CHECK_EQUAL(it->is_valid(), true);
    ++index;
  }
  //std::cout << "collection.size(): " << collection.size() << std::endl;
  BOOST_CHECK_EQUAL(collection.size(), 19);

  /*
  for(auto &dtc : collection) {
    std::cout << "From config" << dtc.sorted().from_config() << std::endl;
    std::cout << "To config " << dtc.sorted().to_config() << std::endl;
  }
  */
}

/*
BOOST_AUTO_TEST_CASE(FCCTest) {

  // Make test project
  BOOST_CHECK_EQUAL(true, true);
  test::FCCTernaryProj proj;
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
  BOOST_CHECK_EQUAL(orbits.size(), 71);

  //print_clust(orbits.begin(), orbits.end(), std::cout, PrototypePrinter<IntegralCluster>());

  // Make PrimPeriodicDiffTransOrbit
  std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
  Kinetics::make_prim_periodic_diff_trans_orbits(
    orbits.begin() + 2,
    orbits.begin() + 4,
    primclex.crystallography_tol(),
    std::back_inserter(diff_trans_orbits),
    &primclex);
  BOOST_CHECK_EQUAL(diff_trans_orbits.size(), 12);

  print_clust(
    diff_trans_orbits.begin(),
    diff_trans_orbits.end(),
    std::cout,
    PrototypePrinter<Kinetics::DiffusionTransformation>());

  Kinetics::DiffusionTransformation diff_trans_prototype = diff_trans_orbits[4].prototype();
  Eigen::Vector3d a1, b1, c1;
  std::tie(a1, b1, c1) = primclex.prim().lattice().vectors();
  Supercell scel {&primclex, Lattice(2 * a1, 2 * b1, 2 * c1)};
  Configuration l12config(scel);
  l12config.init_occupation();
  BOOST_CHECK_EQUAL(true, true);

  //std::cout << l12config << std::endl;
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

  /// Make local clusters
  fs::path l12_local_bspecs_path = "tests/unit/kinetics/l12_local_bspecs_0.json";
  jsonParser l12_local_bspecs {l12_local_bspecs_path};
  std::vector<LocalIntegralClusterOrbit> local_orbits;
  make_local_orbits(
    diff_trans_orbits[4].prototype(),
    l12_local_bspecs,
    alloy_sites_filter,
    primclex.crystallography_tol(),
    std::back_inserter(local_orbits),
    std::cout);
  BOOST_CHECK_EQUAL(true, true);

  /// Check for neighborhood overlap
  BOOST_CHECK_EQUAL(
    Kinetics::has_local_bubble_overlap(local_orbits, l12config.supercell()),
    true);

  /// Constructor enumerator
  Kinetics::DiffTransConfigEnumOccPerturbations enumerator(l12config, diff_trans_orbits[4], l12_local_bspecs);
  BOOST_CHECK_EQUAL(true, true);

  /// Enumerate perturbations (may be duplicates at this point)
  std::vector<Kinetics::DiffTransConfiguration> collection {
    enumerator.begin(),
    enumerator.end()};
  BOOST_CHECK_EQUAL(collection.size(), 1);

  for(auto &dtc : collection) {
    std::cout << "From config" << dtc.sorted().from_config() << std::endl;
    std::cout << "To config " << dtc.sorted().to_config() << std::endl;
  }

}
*/

BOOST_AUTO_TEST_SUITE_END();
