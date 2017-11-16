#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/symmetry/InvariantSubgroup.hh"
#include "casm/symmetry/InvariantSubgroup_impl.hh"
#include "casm/symmetry/Orbit.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/symmetry/SymBasisPermute.hh"

/// What is being used to test it:
#include "Common.hh"
#include "TestConfiguration.hh"
#include "casm/casm_io/jsonFile.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/kinetics/DiffusionTransformation_impl.hh"
#include "casm/symmetry/ConfigSubOrbits_impl.hh"

using namespace CASM;

namespace {
  struct TestConfig0 : test::TestConfiguration {

    TestConfig0(const PrimClex &primclex) :
      TestConfiguration(
        primclex,
        Eigen::Vector3i(2, 1, 1).asDiagonal(), {
      0, 0,  0, 0,  1, 1,  0, 0
    }) {

      BOOST_CHECK_EQUAL(this->scel_fg().size(), 16);
      BOOST_CHECK_EQUAL(this->config_sym_fg().size(), 8);
    }

  };
}

BOOST_AUTO_TEST_SUITE(OrbitTest)

BOOST_AUTO_TEST_CASE(Test0) {
  test::ZrOProj proj;
  proj.check_init();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);
  const Structure &prim = primclex.prim();

  const auto &g = prim.factor_group();

  {
    IntegralCluster generating_element(prim);
    generating_element.elements().push_back(UnitCellCoord(prim, 0, 0, 0, 0));
    PrimPeriodicSymCompare<IntegralCluster> sym_compare(primclex);
    PrimPeriodicOrbit<IntegralCluster> orbit(generating_element, g, sym_compare);
    BOOST_CHECK_EQUAL(orbit.size(), 2);
  }

  {
    IntegralCluster generating_element(prim);
    generating_element.elements().push_back(UnitCellCoord(prim, 1, 0, 0, 0));
    PrimPeriodicSymCompare<IntegralCluster> sym_compare(primclex);
    PrimPeriodicOrbit<IntegralCluster> orbit(generating_element, g, sym_compare);
    BOOST_CHECK_EQUAL(orbit.size(), 2);
  }

  {
    IntegralCluster generating_element(prim);
    generating_element.elements().push_back(UnitCellCoord(prim, 2, 0, 0, 0));
    PrimPeriodicSymCompare<IntegralCluster> sym_compare(primclex);
    PrimPeriodicOrbit<IntegralCluster> orbit(generating_element, g, sym_compare);
    BOOST_CHECK_EQUAL(orbit.size(), 2);
  }

  {
    IntegralCluster generating_element(prim);
    generating_element.elements().push_back(UnitCellCoord(prim, 2, 0, 0, 0));
    generating_element.elements().push_back(UnitCellCoord(prim, 3, 0, 0, 0));
    PrimPeriodicSymCompare<IntegralCluster> sym_compare(primclex);
    PrimPeriodicOrbit<IntegralCluster> orbit(generating_element, g, sym_compare);
    BOOST_CHECK_EQUAL(orbit.size(), 2);
  }

  {
    IntegralCluster generating_element(prim);
    generating_element.elements().push_back(UnitCellCoord(prim, 0, 0, 0, 0));
    generating_element.elements().push_back(UnitCellCoord(prim, 1, 0, 1, 0));
    PrimPeriodicSymCompare<IntegralCluster> sym_compare(primclex);
    PrimPeriodicOrbit<IntegralCluster> orbit(generating_element, g, sym_compare);
    BOOST_CHECK_EQUAL(orbit.size(), 6);
  }


  // --- DiffTrans tests ---

  {
    using namespace Kinetics;

    // Construct
    DiffusionTransformation diff_trans(prim);
    BOOST_CHECK_EQUAL(true, true);
    BOOST_CHECK_EQUAL(diff_trans.occ_transform().size(), 0);

    UnitCellCoord uccoordA(prim, 2, 0, 0, 0);
    UnitCellCoord uccoordB(prim, 3, 0, 0, 0);
    Index iVa = 0;
    Index iO = 1;

    // Add transform (so that it's not sorted as constructed)
    diff_trans.occ_transform().emplace_back(uccoordB, iO, iVa);
    diff_trans.occ_transform().emplace_back(uccoordA, iVa, iO);
    BOOST_CHECK_EQUAL(true, true);
    BOOST_CHECK_EQUAL(diff_trans.is_valid_occ_transform(), true);

    // Add trajectory
    diff_trans.specie_traj().emplace_back(SpecieLocation(uccoordA, iVa, 0), SpecieLocation(uccoordB, iVa, 0));
    diff_trans.specie_traj().emplace_back(SpecieLocation(uccoordB, iO, 0), SpecieLocation(uccoordA, iO, 0));
    BOOST_CHECK_EQUAL(true, true);
    BOOST_CHECK_EQUAL(diff_trans.is_valid_specie_traj(), true);
    BOOST_CHECK_EQUAL(diff_trans.is_valid(), true);

    PrimPeriodicDiffTransSymCompare sym_compare(primclex);
    PrimPeriodicDiffTransOrbit orbit(diff_trans, g, sym_compare, &primclex);
    BOOST_CHECK_EQUAL(orbit.size(), 2);
  }
}

BOOST_AUTO_TEST_CASE(Test1) {

  test::ZrOProj proj;
  proj.check_init();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);
  const Structure &prim = primclex.prim();
  const Lattice &lat = prim.lattice();
  Supercell prim_scel(&primclex, Eigen::Matrix3i::Identity());
  BOOST_CHECK_EQUAL(true, true);

  // Make PrimPeriodicIntegralClusterOrbit
  jsonFile bspecs {"tests/unit/kinetics/ZrO_bspecs_0.json"};

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

  // Make cluster groups & check size, based on reduced symmetry of a vol 2 Supercell
  {
    // - these values have not been checked for correctness, they just check for consistency
    std::vector<Index> orbit_size;
    std::vector<Index> expected_orbit_size = {1, 2, 2, 6, 12, 2, 6, 12, 12, 6, 12,
                                              6, 6, 2, 12, 12, 12, 24, 4, 24, 12, 12, 2, 12, 12, 24, 4, 24, 24, 24, 12,
                                              12, 24, 24, 24, 12, 12, 24, 24, 12, 24, 6, 12, 24, 12, 24, 24, 4, 24, 24,
                                              24, 24, 24, 24, 24, 24, 24, 12, 12, 12, 6, 12, 12, 12, 12, 12, 12, 6, 24,
                                              24, 12, 12, 4, 8
                                             };

    // - these values have not been checked for correctness, they just check for consistency
    std::vector<Index> num_suborbits;
    std::vector<Index> expected_num_suborbits = {1, 2, 1, 4, 3, 2, 4, 6, 4, 4, 3,
                                                 4, 4, 1, 6, 3, 6, 6, 2, 6, 4, 4, 2, 4, 4, 6, 2, 6, 6, 6, 4, 4, 6, 6, 6, 4,
                                                 4, 6, 6, 4, 6, 4, 4, 6, 4, 6, 6, 2, 6, 6, 6, 6, 6, 6, 6, 6, 6, 4, 4, 4, 2,
                                                 4, 4, 3, 3, 4, 4, 4, 6, 6, 4, 4, 2, 2
                                                };

    // - these values have not been checked for correctness, they just check for consistency
    std::vector<Index> suborbit_size;
    std::vector<Index> expected_suborbit_size = {1, 2, 2, 4, 2, 4, 4, 2, 8, 8, 8,
                                                 2, 2, 2, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 8, 8, 4, 2, 4, 4, 2, 8, 8, 8, 2, 4,
                                                 4, 2, 2, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 4, 4, 4, 4, 4, 4, 8, 8, 8,
                                                 8, 8, 8, 4, 4, 8, 8, 8, 8, 8, 8, 4, 8, 8, 4, 4, 8, 8, 4, 2, 2, 4, 8, 8, 4,
                                                 4, 8, 8, 4, 8, 8, 8, 8, 8, 8, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                                 8, 8, 8, 8, 8, 4, 8, 8, 4, 4, 8, 8, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                                 8, 8, 8, 8, 8, 8, 4, 8, 8, 4, 4, 8, 8, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                                 8, 4, 8, 8, 4, 8, 8, 8, 8, 8, 8, 2, 4, 4, 2, 4, 8, 8, 4, 8, 8, 8, 8, 8, 8,
                                                 4, 8, 8, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 4, 4, 8, 8, 8, 8, 8, 8, 8,
                                                 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                                 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 4, 8, 8,
                                                 4, 4, 8, 8, 4, 4, 8, 8, 4, 8, 4, 4, 8, 8, 4, 4, 8, 8, 4, 8, 8, 8, 8, 8, 8,
                                                 4, 8, 8, 4, 4, 8, 8, 4, 2, 4, 4, 2, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 4,
                                                 8, 8, 4, 4, 8, 8, 4, 4, 4, 8, 8
                                                };

    // - these values have not been checked for correctness, they just check for consistency
    std::vector<Index> cluster_group_size;
    std::vector<Index> expected_cluster_group_size = {24, 12, 12, 4, 2, 12, 4,
                                                      2, 2, 4, 2, 4, 4, 12, 2, 2, 2, 1, 6, 1, 2, 2, 12, 2, 2, 1, 6, 1, 1, 1, 2,
                                                      2, 1, 1, 1, 2, 2, 1, 1, 2, 1, 4, 2, 1, 2, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1,
                                                      1, 2, 2, 2, 4, 2, 2, 2, 2, 2, 2, 4, 1, 1, 2, 2, 6, 3
                                                     };

    // Make vol 2 supercell & and background configuration test data ('td')
    TestConfig0 td(primclex);

    Index index = 0;
    Index suborbit_index = 0;
    for(const auto &orbit : orbits) {

      //      std::cout << "\n ---------- Orbit: " << index << " ----------- \n" << std::endl;
      //      std::cout << "orbit.prototype(): \n" << orbit.prototype() << std::endl;
      //      std::cout << "orbit.size(): " << orbit.size() << std::endl;

      orbit_size.push_back(orbit.size());
      BOOST_CHECK_EQUAL(expected_orbit_size[index], orbit.size());

      {
        SymGroup cluster_group = make_invariant_subgroup(orbit);

        //        std::cout << "  cluster_group.size(): " << cluster_group.size() << std::endl;
        cluster_group_size.push_back(cluster_group.size());
        BOOST_CHECK_EQUAL(expected_cluster_group_size[index], cluster_group.size());
      }

      // Test make_invariant_subgroup using orbit generators
      {
        BOOST_CHECK_EQUAL(true, true);
        std::vector<IntegralCluster> generators;
        make_suborbit_generators(orbit, td.config, std::back_inserter(generators));
        BOOST_CHECK_EQUAL(true, true);

        //        std::cout << "suborbits: " << generators.size() << std::endl;
        num_suborbits.push_back(generators.size());
        BOOST_CHECK_EQUAL(expected_num_suborbits[index], generators.size());

        Index N_clust = 0;
        for(const auto &el : generators) {

          BOOST_CHECK_EQUAL(true, true);
          ScelPeriodicIntegralClusterOrbit scel_suborbit(el, td.config_sym_fg(), td.scel_sym_compare);
          BOOST_CHECK_EQUAL(true, true);
          N_clust += scel_suborbit.size();

          BOOST_CHECK_EQUAL(true, true);
          SymGroup cluster_group = make_invariant_subgroup(scel_suborbit);

          BOOST_CHECK_EQUAL(true, true);
          SymGroup cluster_group_alt1 = make_invariant_subgroup(
                                          el, td.config_sym_fg(), td.scel_sym_compare);

          BOOST_CHECK_EQUAL(true, true);
          std::vector<PermuteIterator> cluster_group_alt2 = make_invariant_subgroup(
                                                              el, td.scel, td.config_permute_fg().begin(), td.config_permute_fg().end());

          BOOST_CHECK_EQUAL(true, true);
          //          std::cout << "suborbit: "
          //            << "  scel_fg().size(): " << td.scel.factor_group().size()*td.scel.volume()
          //            << "  td.config_sym_fg().size(): " << td.config_sym_fg().size()
          //            << "  cluster_group.size(): " << cluster_group.size()
          //            << "  cluster_group_alt1.size(): " << cluster_group_alt1.size()
          //            << "  cluster_group_alt2.size(): " << cluster_group_alt2.size()
          //            << "  scel_suborbit.size(): " << scel_suborbit.size()
          //            << std::endl;

          suborbit_size.push_back(scel_suborbit.size());
          BOOST_CHECK_EQUAL(expected_suborbit_size[suborbit_index], scel_suborbit.size());

          BOOST_CHECK_EQUAL(cluster_group.size()*scel_suborbit.size(), td.config_sym_fg().size());
          BOOST_CHECK_EQUAL(cluster_group_alt1.size()*scel_suborbit.size(), td.config_sym_fg().size());
          BOOST_CHECK_EQUAL(cluster_group_alt2.size()*scel_suborbit.size(), td.config_sym_fg().size());

          suborbit_index++;
        }

        if(orbit.prototype().size()) {
          // except for null cluster, check that number of subcluster elements
          //   is consistent with number of orbit elements
          BOOST_CHECK_EQUAL(N_clust, orbit.size()*td.scel.volume());
        }
      }

      index++;
    }

    //    test::print_computed_result(std::cout, "orbit_size", orbit_size);
    //    test::print_computed_result(std::cout, "cluster_group_size", cluster_group_size);
    //    test::print_computed_result(std::cout, "suborbit_size", suborbit_size);
    //    test::print_computed_result(std::cout, "num_suborbits", num_suborbits);
    BOOST_CHECK_EQUAL(true, true);
  }
}

BOOST_AUTO_TEST_SUITE_END()
