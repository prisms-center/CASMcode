#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/clusterography/IntegralCluster.hh"

/// What is being used to test it:
#include "casm/crystallography/Structure.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/app/AppIO.hh"
#include "Common.hh"

using namespace CASM;
using namespace test;

template<typename OrbitType>
jsonParser expected_first_prim_periodic_orbit(const OrbitType &orbit) {
  jsonParser json = jsonParser::object();
  json["equiv"] = jsonParser::array();
  for(const auto &equiv : orbit) {
    json["equiv"].push_back(equiv[0]);
  }
  return json;
}

template<typename OrbitIterator>
jsonParser expected_asym_unit(OrbitIterator begin, OrbitIterator end) {
  jsonParser json = jsonParser::object();
  json["orbits"] = jsonParser::array(std::distance(begin, end), jsonParser::object());
  auto j_it = json["orbits"].begin();
  for(auto it = begin; it != end; ++it) {
    const auto &proto = it->prototype();
    auto &j = (*j_it++)["prototype"];
    j["min_length"] = proto.min_length();
    j["max_length"] = proto.max_length();
    j["sites"].put_array(proto.begin(), proto.end());
  }
  return json;
}

template<typename OrbitIterator>
jsonParser expected_Nclusters(OrbitIterator begin, OrbitIterator end) {
  jsonParser nclust = jsonParser::array();
  // check Nclusters for each orbit
  for(auto it = begin; it != end; ++it) {
    nclust.push_back((int)it->size());
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
    double crystallography_tol = TOL;

    // generate a one site orbit, prim periodic
    {
      IntegralCluster clust(prim);
      BOOST_CHECK_MESSAGE(true, "IntegralCluster constructed");
      clust.elements().push_back(UnitCellCoord(prim, 0, UnitCell(0, 0, 0)));
      BOOST_CHECK_MESSAGE(clust.size() == 1, "site added");
      PrimPeriodicIntegralClusterOrbit orbit(
        clust,
        prim.factor_group(),
        PrimPeriodicIntegralClusterSymCompare(crystallography_tol));
      BOOST_CHECK_MESSAGE(orbit.prototype().size() == 1, "orbit generated");

      check("first_prim_periodic_orbit", j, expected_first_prim_periodic_orbit(orbit), test_cases_path, quiet);
    }

    // generate asym unit
    {
      std::vector<PrimPeriodicIntegralClusterOrbit> asym_unit;
      make_prim_periodic_asymmetric_unit(
        prim,
        alloy_sites_filter,
        crystallography_tol,
        std::back_inserter(asym_unit),
        std::cout);

      // run checks:
      check("asym_unit", j, expected_asym_unit(asym_unit.begin(), asym_unit.end()), test_cases_path, quiet);
    }

    // generate cluster orbits
    {
      std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
      make_prim_periodic_orbits(
        prim,
        j["bspecs"],
        alloy_sites_filter,
        crystallography_tol,
        std::back_inserter(orbits),
        std::cout);

      // run checks:
      check("Nclusters", j, expected_Nclusters(orbits.begin(), orbits.end()), test_cases_path, quiet);
    }

    // ... add more here ...

  }

}

BOOST_AUTO_TEST_SUITE_END()
