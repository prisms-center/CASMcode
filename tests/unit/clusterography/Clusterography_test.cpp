#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/clusterography/IntegralCluster.hh"

/// What is being used to test it:
#include "casm/crystallography/Structure.hh"
#include "casm/clex/PrimClex.hh"
#include "Common.hh"

using namespace CASM;
using namespace test;

template<typename OrbitIterator>
jsonParser expected_Nclusters(OrbitIterator begin, OrbitIterator end) {
  jsonParser nclust = jsonParser::array();
  // check Nclusters for each orbit
  for(auto it = begin; it != end; ++it) {
    nclust.push_back(it->size());
  }
  return nclust;
}

BOOST_AUTO_TEST_SUITE(BasicStructureSiteTest)

BOOST_AUTO_TEST_CASE(ClusterographyTest) {

  // read test file
  fs::path test_cases_path("tests/unit/clusterography/test_cases.json");
  jsonParser tests(test_cases_path);

  for(auto test_it = tests.begin(); test_it != tests.end(); ++test_it) {

    // input and expected output data
    jsonParser &j = *test_it;

    // if false: print calculated results if no test data; if true: suppress
    bool quiet = false;
    j.get_else(quiet, "quiet", false);

    BOOST_CHECK_MESSAGE(j.contains("title"), "test case 'title' is required");
    BOOST_CHECK_MESSAGE(j.contains("prim"), "test case 'prim' is required");
    BOOST_CHECK_MESSAGE(j.contains("bspecs"), "test case 'bspecs' is required");

    // generate prim
    Structure prim(read_prim(j["prim"]));

    // generate cluster orbits
    std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
    make_orbits(prim,
                prim.factor_group(),
                j["bspecs"],
                TOL,
                alloy_sites_filter,
                PrimPeriodicIntegralClusterSymCompare(TOL),
                std::back_inserter(orbits),
                std::cout);

    // run checks:
    check("Nclusters", j, expected_Nclusters(orbits.begin(), orbits.end()), test_cases_path, quiet);
    // ... add more here ...

  }

}

BOOST_AUTO_TEST_SUITE_END()
