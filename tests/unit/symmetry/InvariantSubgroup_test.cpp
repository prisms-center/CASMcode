#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/symmetry/InvariantSubgroup.hh"
#include "casm/symmetry/InvariantSubgroup_impl.hh"

/// What is being used to test it:
#include "Common.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/clusterography/ClusterOrbits.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(InvariantSubgroupTest)

BOOST_AUTO_TEST_CASE(Test0) {
  test::ZrOProj proj;
  proj.check_init();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);
  const Structure &prim = primclex.prim();
  const Lattice &lat = prim.lattice();
  Supercell prim_scel(&primclex, Eigen::Matrix3i::Identity());

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

  // Make cluster groups & check size, based on prim.factor_group symmetry
  {
    // - these values have not been checked for correctness, they just check for consistency
    std::vector<Index> expected_cluster_group_size = {24, 12, 12, 4, 2, 12, 4, 2,
                                                      2, 4, 2, 4, 4, 12, 2, 2, 2, 1, 6, 1, 2, 2, 12, 2, 2, 1, 6, 1, 1, 1, 2, 2, 1,
                                                      1, 1, 2, 2, 1, 1, 2, 1, 4, 2, 1, 2, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
                                                      2, 2, 4, 2, 2, 2, 2, 2, 2, 4, 1, 1, 2, 2, 6, 3
                                                     };
    std::vector<Index> cluster_group_size;
    Index index = 0;
    for(const auto &orbit : orbits) {

      // Test make_invariant_subgroup using orbit generators
      {
        SymGroup cluster_group = make_invariant_subgroup(
                                   orbit.prototype(),
                                   primclex.prim().factor_group(),
                                   orbit.sym_compare());
        cluster_group_size.push_back(cluster_group.size());
        BOOST_CHECK_EQUAL(expected_cluster_group_size[index], cluster_group.size());
      }

      // Test make_invariant_subgroup using the orbit equivalence map & orbit.prototype()
      {
        SymGroup cluster_group = make_invariant_subgroup(orbit);
        BOOST_CHECK_EQUAL(expected_cluster_group_size[index], cluster_group.size());
      }

      // Test make_invariant_subgroup using Supercell
      {
        std::vector<PermuteIterator> cluster_group = make_invariant_subgroup(orbit.prototype(), prim_scel);
        BOOST_CHECK_EQUAL(expected_cluster_group_size[index], cluster_group.size());
      }

      index++;
    }
    //test::print_computed_result(std::cout, "cluster_group_size", cluster_group_size);
    BOOST_CHECK_EQUAL(true, true);
  }

  // Make cluster groups & check size, based on reduced symmetry of a vol 2 Supercell
  {
    // Make vol 2 supercell & and background configuration
    Eigen::Vector3d a, b, c;
    std::tie(a, b, c) = primclex.prim().lattice().vectors();
    Supercell scel_vol2 {&primclex, Lattice(2 * a, 1 * b, 1 * c)};
    Configuration config(scel_vol2);

    // Scel sym_compare
    ScelPeriodicSymCompare<IntegralCluster> scel_sym_compare(
      scel_vol2.prim_grid(),
      scel_vol2.crystallography_tol());


    // Get the config factor group (should just be all Supercell operations)
    std::vector<PermuteIterator> _config_fg = config.factor_group();
    SymGroup config_fg = make_sym_group(_config_fg.begin(), _config_fg.end());

    // - these values have not been checked for correctness, they just check for consistency
    std::vector<Index> expected_cluster_group_size = {16, 4, 4, 4, 2, 4, 2, 2, 1,
                                                      4, 2, 2, 2, 4, 2, 2, 2, 1, 2, 1, 1, 2, 4, 2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1,
                                                      1, 1, 2, 1, 1, 1, 1, 4, 2, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1,
                                                      2, 4, 2, 1, 2, 2, 2, 1, 2, 1, 1, 1, 2, 2, 1
                                                     };
    std::vector<Index> cluster_group_size;
    Index index = 0;
    for(const auto &orbit : orbits) {

      // Test make_invariant_subgroup using orbit generators
      {
        SymGroup cluster_group = make_invariant_subgroup(
                                   orbit.prototype(),
                                   config_fg,
                                   scel_sym_compare);
        cluster_group_size.push_back(cluster_group.size());
        BOOST_CHECK_EQUAL(expected_cluster_group_size[index], cluster_group.size());
      }

      // Test make_invariant_subgroup using the orbit equivalence map & orbit.prototype()
      {
        Orbit<IntegralCluster, ScelPeriodicSymCompare<IntegralCluster>> suborbit(
                                                                       orbit.prototype(), config_fg, scel_sym_compare);
        SymGroup cluster_group = make_invariant_subgroup(suborbit);
        BOOST_CHECK_EQUAL(expected_cluster_group_size[index], cluster_group.size());
      }

      // Test make_invariant_subgroup using Supercell
      {
        std::vector<PermuteIterator> cluster_group = make_invariant_subgroup(orbit.prototype(), scel_vol2);
        BOOST_CHECK_EQUAL(expected_cluster_group_size[index], cluster_group.size());
      }

      index++;
    }
    //test::print_computed_result(std::cout, "cluster_group_size", cluster_group_size);
    BOOST_CHECK_EQUAL(true, true);
  }
}

BOOST_AUTO_TEST_SUITE_END()
