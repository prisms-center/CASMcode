#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clusterography/Orbitree.hh"

/// What is being used to test it:
#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/app/AppIO.hh"
#include "casm/clusterography/jsonClust.hh"
#include "Common.hh"

using namespace CASM;
using namespace test;

jsonParser expected_Nclusters(const SiteOrbitree &tree) {
  jsonParser nclust = jsonParser::array();
  // check Nclusters for each orbit
  for(int branch = 0; branch < tree.size(); ++branch) {
    nclust.push_back(jsonParser::array());
    for(int orbit = 0; orbit < tree[branch].size(); ++orbit) {
      nclust[branch].push_back(tree[branch][orbit].size());
    }
  }
  return nclust;
}

BOOST_AUTO_TEST_SUITE(BasicStructureSiteTest)

BOOST_AUTO_TEST_CASE(ClusterographyTest) {

  // read test file
  fs::path test_cases_path("tests/unit/clusterography/test_cases.json");
  jsonParser tests(test_cases_path);
  double tol = TOL;

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
    SiteOrbitree tree(prim.lattice(), tol);
    tree = make_orbitree(prim, j["bspecs"], tol);

    jsonParser clust_json;
    to_json(jsonHelper(tree, prim), clust_json);

    // run checks:
    check("Nclusters", j, expected_Nclusters(tree), test_cases_path, quiet);
    // ... add more here ...

  }

}

BOOST_AUTO_TEST_SUITE_END()
