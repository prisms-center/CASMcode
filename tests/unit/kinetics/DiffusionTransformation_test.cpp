#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffusionTransformationEnum.hh"
#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"

/// What is being used to test it:
#include "casm/clex/PrimClex.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/AppIO_impl.hh"
#include "Common.hh"
#include "casm/completer/Handlers.hh"
#include "casm/app/casm_functions.hh"
#include "casm/app/enum.hh"

using namespace CASM;
using namespace test;

BOOST_AUTO_TEST_SUITE(DiffusionTransformationTest)

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

  Kinetics::DiffusionTransformation trans(primclex.prim());
  BOOST_CHECK_EQUAL(trans.size(), 0);

  // in ZrO (and other binary Va-X materials):
  // empty and point cluster have 0 valid DiffusionTransformation
  // 2-cluster have 3 valid (Va-O, O-Va, O-O)
  // 3-cluster have 8 valid (Va-O-O & rev, O-Va-O & rev, O-O-Va & rev, O-O-O & rev)
  // 4-cluster have 69 valid
  //   4 (Va-O-O-O) * 9 (6 complete perm, 3 simultaneous pair perm) = 36
  //   6 (Va-Va-O-O) * 4 (2 complete perm, 2 simultaneous pair perm) = 24
  //   4 (Va-Va-Va-O) * 0 = 0
  //   1 (O-O-O-O) * 9 = 9
  std::map<int, int> count_check = {{0, 0}, {1, 0}, {2, 3}, {3, 8}, {4, 69}};

  // these have not all been manually checked
  std::vector<int> orbit_count;
  std::vector<int> expected_orbit_count = {0, 0, 2, 2, 2, 2, 2, 2, 2,
                                           2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 3, 3, 3, 3, 3, 4, 2, 4, 4, 4, 3, 3, 4, 4, 4,
                                           4, 4, 4, 4, 4, 4, 3, 3, 4, 3, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 27,
                                           20, 20, 30, 30, 20, 20, 17, 36, 36, 20, 30, 8, 12
                                          };
  auto expected_orbit_count_it = expected_orbit_count.begin();

  // these have not all been manually checked
  std::vector<int> mult;
  std::vector<int> expected_mult = {2, 2, 6, 6, 12, 12, 2, 2, 6, 6, 12, 12, 12,
                                    12, 6, 6, 12, 12, 6, 6, 6, 6, 2, 2, 12, 12, 12, 12, 12, 12, 24, 24, 12, 4, 24,
                                    24, 24, 24, 12, 24, 12, 24, 12, 12, 2, 4, 2, 24, 12, 12, 12, 24, 12, 24, 24,
                                    24, 24, 12, 4, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 12, 24, 12, 24,
                                    12, 12, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 12, 12, 12, 12, 12,
                                    12, 12, 12, 24, 24, 24, 24, 24, 24, 24, 24, 12, 12, 12, 12, 24, 24, 24, 24, 6,
                                    12, 6, 12, 24, 12, 24, 24, 24, 24, 12, 24, 12, 24, 24, 24, 24, 24, 24, 24, 24,
                                    12, 4, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
                                    24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
                                    12, 12, 12, 24, 12, 24, 12, 12, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 12, 12,
                                    12, 12, 12, 12, 12, 12, 12, 6, 6, 6, 6, 6, 6, 24, 24, 24, 24, 24, 24, 12, 24,
                                    24, 24, 24, 24, 24, 24, 12, 24, 12, 24, 12, 24, 24, 24, 24, 24, 24, 24, 12,
                                    24, 24, 24, 24, 24, 24, 24, 12, 24, 12, 24, 12, 24, 12, 12, 12, 12, 12, 12,
                                    12, 12, 12, 12, 12, 12, 12, 12, 12, 24, 24, 24, 24, 24, 24, 12, 12, 12, 12,
                                    12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
                                    12, 24, 24, 24, 24, 24, 24, 12, 12, 12, 12, 12, 12, 12, 12, 12, 24, 24, 24,
                                    24, 24, 24, 12, 24, 24, 24, 24, 12, 24, 24, 24, 24, 12, 12, 24, 24, 24, 24,
                                    24, 24, 24, 24, 12, 24, 24, 24, 24, 24, 24, 24, 12, 24, 12, 24, 12, 24, 12,
                                    12, 12, 12, 12, 12, 12, 6, 12, 6, 24, 24, 24, 6, 12, 12, 6, 24, 24, 24, 24,
                                    24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
                                    24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
                                    24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
                                    24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 12,
                                    24, 24, 24, 24, 12, 24, 24, 24, 12, 24, 24, 12, 12, 12, 12, 12, 12, 12, 12,
                                    12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 24, 24, 24, 24, 24, 24, 12, 12,
                                    12, 12, 12, 12, 24, 24, 12, 24, 12, 24, 12, 12, 24, 24, 24, 24, 24, 24, 24,
                                    24, 24, 24, 24, 24
                                   };
  auto expected_mult_it = expected_mult.begin();


  for(auto it = orbits.begin(); it != orbits.end(); ++it) {

    // print the IntegralCluster prototype used to generate DiffTrans
    //std::cout << "----------------\n";
    //std::cout << "IntegralCluster: \n" << it->prototype() << std::endl;

    // check the number of valid DiffTrans
    Kinetics::DiffusionTransformationEnum e {it->prototype()};
    BOOST_CHECK_EQUAL(std::distance(e.begin(), e.end()), count_check[it->prototype().size()]);

    // make orbits of DiffTrans
    std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
    Kinetics::make_prim_periodic_diff_trans_orbits(it, it + 1, primclex.crystallography_tol(), std::back_inserter(diff_trans_orbits));

    // check how many DiffTrans orbits there are for each IntegralCluster orbit
    orbit_count.push_back(diff_trans_orbits.size());
    if(expected_orbit_count_it != expected_orbit_count.end()) {
      BOOST_CHECK_EQUAL(diff_trans_orbits.size(), *expected_orbit_count_it);
      ++expected_orbit_count_it;
    }
    else {
      BOOST_CHECK_EQUAL(diff_trans_orbits.size(), 0);
    }

    // check the size of the DiffTrans orbits
    for(const auto &orb : diff_trans_orbits) {
      mult.push_back(orb.size());
      if(expected_mult_it != expected_mult.end()) {
        BOOST_CHECK_EQUAL(orb.size(), *expected_mult_it);
        ++expected_mult_it;
      }
      else {
        BOOST_CHECK_EQUAL(orb.size(), 0);
      }
    }
    jsonParser dtjson;
    dtjson.put_obj();
    if(diff_trans_orbits.size()) {
      to_json(diff_trans_orbits[0].prototype(), dtjson);
      Kinetics::DiffusionTransformation trans = jsonConstructor<Kinetics::DiffusionTransformation>::from_json(dtjson, primclex.prim());
      BOOST_CHECK_EQUAL(trans == diff_trans_orbits[0].prototype(), 1);
      BOOST_CHECK_EQUAL(diff_trans_orbits[0].prototype().is_valid(), 1);
      BOOST_CHECK_EQUAL(trans.is_valid(), 1);
    }

    //print DiffTrans prototypes
    //{
    //  PrototypePrinter<Kinetics::DiffusionTransformation> printer;
    //  print_clust(diff_trans_orbits.begin(), diff_trans_orbits.end(), std::cout, printer);
    //}
  }

  fs::path difftrans_path = "tests/unit/kinetics/diff_trans.json";
  jsonParser diff_trans_json {difftrans_path};
  Completer::EnumOption enum_opt;
  enum_opt.desc();
  int success = Kinetics::DiffusionTransformationEnum::run(primclex, diff_trans_json, enum_opt);
  BOOST_CHECK_EQUAL(success, 0);

  /*
  auto vecprinter = [=](const std::string name, const std::vector<int>& v) {
    std::cout << name << " = {";
    for(int i=0; i<v.size(); ++i) {
      std::cout << v[i];
      if(i != v.size()-1) {
        std::cout << ", ";
      }
    }
    std::cout << "};" << std::endl;
  };

  vecprinter("orbit_count", orbit_count);
  vecprinter("mult", mult);
  */

}

BOOST_AUTO_TEST_SUITE_END()
