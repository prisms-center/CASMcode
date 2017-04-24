#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/kinetics/DiffTransConfigEnumPerturbations.hh"
#include "casm/kinetics/DiffTransEnumEquivalents.hh"

/// What is being used to test it:
#include "casm/clex/PrimClex.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/AppIO_impl.hh"
#include "Common.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/kinetics/DoFTransformation.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffusionTransformationEnum.hh"
#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"
#include "casm/clusterography/ClusterOrbits.hh"

using namespace CASM;
using namespace test;

BOOST_AUTO_TEST_SUITE(DiffTransConfigEnumPerturbationsTest)

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


  /// Make background_config
  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();
  Supercell scel {&primclex, Lattice(2 * a, 2 * b, 3 * c)};
  Configuration config(scel);
  config.init_occupation();
  config.init_displacement();
  config.init_deformation();
  config.init_specie_id();
  config.set_occupation({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0});
  //  config.print_occupation(std::cout);
  //  std::cout << config;
  Configuration bg_config_prim = config.primitive();

  //  std::cout << "Supercell factor groups" << "\n";
  //  for ( auto &g: bg_config_prim.supercell().factor_group()){
  //  g.print(std::cout);
  //  }

  /// Find prototype of m_diff_trans_orbit
  //print_clust(orbits.begin() + 2, orbits.begin() + 3, std::cout, ProtoSitesPrinter());
  std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
  Kinetics::make_prim_periodic_diff_trans_orbits(orbits.begin() + 2, orbits.begin() + 3, primclex.crystallography_tol(), std::back_inserter(diff_trans_orbits));
  Kinetics::DiffusionTransformation diff_trans_prototype = diff_trans_orbits[0].prototype();
  std::cout << diff_trans_orbits.size() << "\n";

  std::cout << "Prototype Diff Trans:" << "\n" << diff_trans_prototype << "\n";

  /// Find unique DiffusionTransformations
//  PermuteIterator begin = bg_config_prim.supercell().permute_begin();
//  PermuteIterator end = bg_config_prim.supercell().permute_end();
//  Kinetics::DiffTransEnumEquivalents diff_trans_unique(diff_trans_prototype, begin, end, bg_config_prim);
  PermuteIterator begin = config.supercell().permute_begin();
  PermuteIterator end = config.supercell().permute_end();
  Kinetics::DiffTransEnumEquivalents diff_trans_unique(diff_trans_prototype, begin, end, config);

  //  std::cout << "Config factor groups" << "\n";
  //  for ( auto &g: diff_trans_unique.invariant_subgroup()){
  //  g.sym_op().print(std::cout);
  //  }

  std::cout << "Diff_Trans_Unique_Current" << "\n";
  std::vector<Kinetics::DiffusionTransformation> subdifftrans;


  for(auto it = diff_trans_unique.begin(); it != diff_trans_unique.end(); ++it) {
    subdifftrans.push_back(*it);
  }

  for(auto it = subdifftrans.begin(); it != subdifftrans.end(); ++it) {
    std::cout << *it << "\n";
  }

  /// Check if number of unique diff trans is what we expect
  //Size of diff trans factor group should be equal to (size of scel fg)*(volume of scel)/(size of bg_config fg)
  std::cout << "size of Config factor group " << diff_trans_unique.invariant_subgroup().size() << std::endl;
  std::cout << "size of Supercell factor group " << bg_config_prim.supercell().factor_group().size() << std::endl;
  std::cout << "size of Supercell (#prims) " << bg_config_prim.supercell().volume() << std::endl;
  std::cout << "size of non-prim Supercell factor group " << config.supercell().factor_group().size() << std::endl;
  std::cout << "size of non-prim Supercell (#prims) " << config.supercell().volume() << std::endl;
  std::cout << "size of diff trans unique " << subdifftrans.size() << std::endl;

  int config_fg_int = static_cast<int>(diff_trans_unique.invariant_subgroup().size());
  int scel_fg_int = static_cast<int>(bg_config_prim.supercell().factor_group().size());
  int scel_vol_int = static_cast<int>(bg_config_prim.supercell().volume());

  int subdifftrans_size_int = static_cast<int>(subdifftrans.size());

  int expected_subdifftrans_size_int = scel_fg_int * scel_vol_int / config_fg_int;
  int expected_subdifftrans_size_int_remainder = scel_fg_int * scel_vol_int % config_fg_int;

  BOOST_CHECK_EQUAL(expected_subdifftrans_size_int_remainder, 0);
  BOOST_CHECK_EQUAL(expected_subdifftrans_size_int, subdifftrans_size_int);

  /// Check that the unique diff trans found are the ones we expect
  // Sites and clusters are equivalent along b lattice vector direction
  // Side refers exclusively to the a lattice vector direction
  // Above is in the positive c lattice direction
  // Below is in the negative c lattice direction
  // Expected unique diff trans indicating pairs by linear index and description:

  //  1
  //  Va: 2, 0 0 0 : 0 0  ->  3, 0 0 0 : 0 0
  //  O: 3, 0 0 0 : 1 0  ->  2, 0 0 0 : 1 0
  //  24-36, 25-37
  //  Ox-Ox, 1 Ox above, 2 Ox below, 1 Va to either side of upper site
  //
  //  2
  //  Va: 2, 1 1 1 : 0 0  ->  3, 1 1 1 : 0 0
  //  O: 3, 1 1 1 : 1 0  ->  2, 1 1 1 : 1 0
  //  27-39, 26-38
  //  Ox-Ox, 1 Ox above, 1 Va below, 1 Va to either side of upper site
  //
  //  3
  //  Va: 2, 0 0 2 : 0 0  ->  3, 0 0 2 : 0 0
  //  O: 3, 0 0 2 : 1 0  ->  2, 0 0 2 : 1 0
  //  28-40, 29-41
  //  Ox-Ox, 3 Ox above, 1 Va below, 1 Va to either side of upper site
  //
  //  4
  //  Va: 2, 1 1 0 : 0 0  ->  3, 1 1 0 : 0 0
  //  O: 3, 1 1 0 : 1 0  ->  2, 1 1 0 : 1 0
  //  31-43, 30-42
  //  Ox-Va, 3 Ox above, 1 Va below, no NN Va to the sides
  //
  //  5
  //  Va: 2, 0 0 1 : 0 0  ->  3, 0 0 1 : 0 0
  //  O: 3, 0 0 1 : 1 0  ->  2, 0 0 1 : 1 0
  //  32-44, 33-45
  //  Ox-Va, 5 Ox above, 4 Ox below, no NN Va to the sides
  //
  //  6
  //  Va: 2, 1 1 2 : 0 0  ->  3, 1 1 2 : 0 0
  //  O: 3, 1 1 2 : 1 0  ->  2, 1 1 2 : 1 0
  //  35-47, 34-46
  //  Ox-Va, 1 Ox above, 2 Ox below, no NN Va to the sides
  //
  //  Clusters 7-12 are similar to 1-6, but are unique due to lack of 2 fold symmetry in config
  //
  //  7
  //  Va: 3,  0 -1  0 : 0 0  ->  2,  0 -1  1 : 0 0
  //  O: 2,  0 -1  1 : 1 0  ->  3,  0 -1  0 : 1 0
  //  37-33, 36-32
  //  Ox-Ox, 1 Va above, 3 Ox below, 1 Va to either side of lower site
  //  Similar to 3
  //
  //  8
  //  Va: 3, 1 0 1 : 0 0  ->  2, 1 0 2 : 0 0
  //  O: 2, 1 0 2 : 1 0  ->  3, 1 0 1 : 1 0
  //  38-34, 39-35
  //  Ox-Ox, 1 Va above, 1 Ox below, 1 Va to either side of lower site
  //  Similar to 2
  //
  //  9
  //  Va: 3,  0 -1  2 : 0 0  ->  2,  0 -1  3 : 0 0
  //  O: 2,  0 -1  3 : 1 0  ->  3,  0 -1  2 : 1 0
  //  41-25, 40-24
  //  Ox-Ox, 2 Ox above, 1 Ox below, 1 Va to either side of lower site
  //  Similar to 1
  //
  //  10
  //  Va: 3, 1 0 0 : 0 0  ->  2, 1 0 1 : 0 0
  //  O: 2, 1 0 1 : 1 0  ->  3, 1 0 0 : 1 0
  //  42-26, 43-27
  //  Va-Ox, 2 Ox above, 1 Ox below, no NN Va to the sides
  //  Similar to 6
  //
  //  11
  //  Va: 3,  0 -1  1 : 0 0  ->  2,  0 -1  2 : 0 0
  //  O: 2,  0 -1  2 : 1 0  ->  3,  0 -1  1 : 1 0
  //  45-29, 44-28
  //  Va-Ox, 4 Ox above, 5 Ox below, no NN Va to the sides
  //  Similar to 5
  //
  //  12
  //  Va: 3, 1 0 2 : 0 0  ->  2, 1 0 3 : 0 0
  //  O: 2, 1 0 3 : 1 0  ->  3, 1 0 2 : 1 0
  //  46-30, 47-31
  //  Va-Ox, 1 Va above, 3 Ox below, no NN Va to the sides
  //  Similar to 4


  for(auto it = subdifftrans.begin(); it != subdifftrans.end(); ++it) {
    //std::vector<Kinetics::SpecieTrajectory> tmp = *it.specie_traj();
    std::vector<Kinetics::SpecieTrajectory> tmp = it->specie_traj();
    for(auto it2 = tmp.begin(); it2 != tmp.end(); ++it2) {
      Kinetics::SpecieLocation tmpfrom = it2->from;
      Kinetics::SpecieLocation tmpto = it2->to;
      UnitCellCoord tmpfromcoord = tmpfrom.uccoord;
      UnitCellCoord tmptocoord = tmpto.uccoord;
      //std::cout << "tmp from coord: " << tmpfromcoord << std::endl;
      //tmpfromcoord.coordinate().print(std::cout);
      //std::cout << "\n";
      //tmpfromcoord.unit().main_print(std::cout, FRAC, true, 0);
    }
  }

  //
  //diff_trans_unique.my_increment();
  //std::cout << "\n";
  //std::cout << diff_trans_unique.current();

  //DiffTransEnumEquivalents diff_trans_enum;
  //diff_trans_enum.current()

}
BOOST_AUTO_TEST_SUITE_END()
