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
  config.print_occupation(std::cout);
  std::cout << config;
  Configuration bg_config_prim = config.primitive();

  /// Find prototype of m_diff_trans_orbit
  //print_clust(orbits.begin() + 2, orbits.begin() + 3, std::cout, ProtoSitesPrinter());
  std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
  Kinetics::make_prim_periodic_diff_trans_orbits(orbits.begin() + 2, orbits.begin() + 3, primclex.crystallography_tol(), std::back_inserter(diff_trans_orbits));
  Kinetics::DiffusionTransformation diff_trans_prototype = diff_trans_orbits[0].prototype();

  std::cout << "Prototype Diff Trans:" << "\n" << diff_trans_prototype << "\n";

  /// Find unique DiffusionTransformations
  PermuteIterator begin = bg_config_prim.supercell().permute_begin();
  PermuteIterator end = bg_config_prim.supercell().permute_end();
  Kinetics::DiffTransEnumEquivalents diff_trans_unique(diff_trans_prototype, begin, end, bg_config_prim);
  std::cout << "size of Config factor group " << diff_trans_unique.invariant_subgroup().size() << std::endl;
  std::cout << "size of Supercell factor group " << bg_config_prim.supercell().factor_group().size() << std::endl;
  std::cout << "size of Supercell (#prims) " << bg_config_prim.supercell().volume() << std::endl;

  std::vector<Kinetics::DiffusionTransformation> subdifftrans;


  for(auto it = diff_trans_unique.begin(); it != diff_trans_unique.end(); ++it) {
    subdifftrans.push_back(*it);
  }
  std::cout << "Diff_Trans_Uniques " << subdifftrans.size() << " total" << "\n";
  for(auto it = subdifftrans.begin(); it != subdifftrans.end(); ++it) {
    std::cout << *it;
    std::cout << "\n";
  }
  //
  //diff_trans_unique.my_increment();
  //std::cout << "\n";
  //std::cout << diff_trans_unique.current();

  //DiffTransEnumEquivalents diff_trans_enum;
  //diff_trans_enum.current()

}
BOOST_AUTO_TEST_SUITE_END()
