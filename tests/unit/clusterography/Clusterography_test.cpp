#include "gtest/gtest.h"
#include "autotools.hh"

#include <memory>

/// What is being tested:
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/clusterography/IntegralCluster.hh"

/// What is being used to test it:
#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/database/Database.hh"
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

TEST(BasicStructureSiteTest, ClusterographyTest) {

  // read test file
  fs::path test_cases_path(autotools::abs_srcdir() + "/tests/unit/clusterography/test_cases.json");
  jsonParser tests(test_cases_path);
  double tol = TOL;

  Log &log = null_log();

  for(auto test_it = tests.begin(); test_it != tests.end(); ++test_it) {

    // input and expected output data
    jsonParser &j = *test_it;

    // if false: print calculated results if no test data; if true: suppress
    bool quiet = false;
    j.get_else(quiet, "quiet", false);

    EXPECT_TRUE(j.contains("title")) << "test case 'title' is required";
    EXPECT_TRUE(j.contains("prim")) << "test case 'prim' is required";
    EXPECT_TRUE(j.contains("bspecs")) << "test case 'bspecs' is required";

    // generate prim
    auto prim_ptr = std::make_shared<Structure>(read_prim(j["prim"], TOL));
    const auto &prim = *prim_ptr;
    double crystallography_tol = TOL;

    // generate a one site orbit, prim periodic
    {
      IntegralCluster clust(prim);
      EXPECT_TRUE(true) << "IntegralCluster constructed";
      clust.elements().emplace_back(0, UnitCell(0, 0, 0));
      EXPECT_TRUE(clust.size() == 1) << "site added";
      PrimPeriodicOrbit<IntegralCluster> orbit(
        clust,
        prim.factor_group(),
        PrimPeriodicSymCompare<IntegralCluster>(prim_ptr, crystallography_tol));
      EXPECT_TRUE(orbit.prototype().size() == 1) << "orbit generated";

      check("first_prim_periodic_orbit", j, expected_first_prim_periodic_orbit(orbit), test_cases_path, quiet);
    }

    // generate asym unit
    {
      std::vector<PrimPeriodicOrbit<IntegralCluster>> asym_unit;
      make_prim_periodic_asymmetric_unit(
        prim,
        alloy_sites_filter,
        crystallography_tol,
        std::back_inserter(asym_unit),
        log);

      // run checks:
      check("asym_unit", j, expected_asym_unit(asym_unit.begin(), asym_unit.end()), test_cases_path, quiet);
    }

    // generate cluster orbits
    {
      std::vector<PrimPeriodicOrbit<IntegralCluster>> orbits;
      make_prim_periodic_orbits(
        prim,
        j["bspecs"],
        alloy_sites_filter,
        crystallography_tol,
        std::back_inserter(orbits),
        log);

      // run checks:
      check("Nclusters", j, expected_Nclusters(orbits.begin(), orbits.end()), test_cases_path, quiet);
    }

    // ... add more here ...

  }

}
