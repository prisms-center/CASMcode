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
#include "casm/kinetics/SubOrbitGenerators.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/casm_io/VaspIO.hh"

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

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();


  /// Make background_config
  Supercell scel {&primclex, Lattice(2 * a, 2 * b, 3 * c)};
  Configuration config(scel);
  config.init_occupation();
  config.init_displacement();
  config.init_deformation();
  config.init_specie_id();
  config.set_occupation({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0});
  //config.set_occupation({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
  //  config.print_occupation(std::cout);
  //  std::cout << config;

  /*
  fs::ofstream file;
  fs::path POSCARpath = "POSCAR";
  file.open(POSCARpath);
  VaspIO::PrintPOSCAR p(config);
  p.sort();
  p.print(file);
  file.close();
  */


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
  //std::cout << diff_trans_orbits.size() << "\n";

  //std::cout << "Prototype Diff Trans:" << "\n" << diff_trans_prototype << "\n";

  /// Test bubble checkers
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

  BOOST_CHECK_EQUAL(Kinetics::has_local_bubble_overlap(local_orbits, scel1), 1);
  BOOST_CHECK_EQUAL(Kinetics::has_local_bubble_overlap(local_orbits, scel2), 0);
  BOOST_CHECK_EQUAL(Kinetics::has_local_bubble_overlap(local_orbits, scel3), 1);
  BOOST_CHECK_EQUAL(Kinetics::has_local_bubble_overlap(local_orbits, scel4), 1);
  std::vector<Supercell> result = Kinetics::viable_supercells(local_orbits, scel_list);
  BOOST_CHECK_EQUAL(*(result.begin()) == scel2, 1);

  //const multivector<SymOp>::X<2> diff_trans_map = diff_trans_orbits[0].equivalence_map();

  /*
  int c1 = 1;
  for(auto it = diff_trans_map.begin(); it != diff_trans_map.end(); ++it) {
    std::cout << c1 << std::endl;
    ++c1;
    for(auto it2 = it->begin(); it2 != it->end(); ++it2) {
      it2->print(std::cout);
    }
  }
  /**/


  //const CASM::SymOp Tk = scel.prim_grid().sym_op(0);
  //Tk.print(std::cout);



  const Orbit<Kinetics::DiffusionTransformation, Kinetics::PrimPeriodicDiffTransSymCompare> orbit = diff_trans_orbits[0];
  const SymGroup config_scel_fg = config.supercell().factor_group();
  const SymGroup prim_config_scel_fg = bg_config_prim.supercell().factor_group();

  /*
  std::vector<SymOp> unique_prim_op;
  for(prim_op : prim_config_scel_fg) {
    if(prim_op <= scel_op*prim_op for all scel_op) {
      unique_prim_op.push_back(prim_op);
    }
  }
  /**/

  const Configuration prim_config = config.primitive();
  std::vector<Kinetics::DiffusionTransformation> unique_prim_config_elements;
  const Supercell prim_scel = prim_config.supercell();
  const std::vector<CASM::PermuteIterator> prim_config_fg = prim_config.factor_group();

  for(auto i = 0; i < orbit.size(); i++) {
    //std::cout << "i: " << i << std::endl;
    for(auto k = 0; k < prim_scel.prim_grid().size(); k++) {
      //std::cout << "k: " << k << std::endl;
      SymOp Sik = prim_scel.prim_grid().sym_op(k) * orbit.equivalence_map()[i][0];
      //std::cout << "Sik for i = " << i << " and k = " << k << "\n";
      //Sik.print(std::cout);
      //std::cout << std::endl;
      auto lambda = [&](const PermuteIterator & it) {
        Kinetics::ScelPeriodicDiffTransSymCompare symcompare(prim_scel.prim_grid(), primclex.crystallography_tol());
        // get primclex from something for template
        return symcompare.compare(copy_apply(Sik, orbit.prototype()), copy_apply(it.sym_op(), copy_apply(Sik, orbit.prototype())));
      };

      if(std::none_of(prim_config_fg.begin(), prim_config_fg.end(), lambda)) {
        Kinetics::DiffusionTransformation Eik = copy_apply(Sik, orbit.prototype());
        unique_prim_config_elements.push_back(Eik);
        //std::cout << Eik << std::endl;
      }
    }
  }

  std::cout << "Number of elements in unique_prim_config_elements: " << unique_prim_config_elements.size() << std::endl;

  // Map the unique orbit elements found above onto config by using SymOps in primcell that are not present in supercell

  std::vector<SymOp> unique_prim_op;

  for(auto prim_op : prim_config_scel_fg) {

    auto lambda = [&](const SymOp & config_op) {
      Kinetics::ScelPeriodicDiffTransSymCompare symcompare(prim_scel.prim_grid(), primclex.crystallography_tol());
      return symcompare.compare(copy_apply(prim_op, orbit.prototype()), copy_apply(config_op, copy_apply(prim_op, orbit.prototype())));
    };

    if(std::none_of(config_scel_fg.begin(), config_scel_fg.end(), lambda)) {
      unique_prim_op.push_back(prim_op);
    }
  }

  std::vector<Kinetics::DiffusionTransformation> unique_config_elements;
  for(auto e : unique_prim_config_elements) {
    for(auto op : unique_prim_op) {
      unique_config_elements.push_back(copy_apply(op, e));
      //*result++ = copy_apply(op, e);
    }
  }

  std::cout << "Number of elements in unique_config_elements: " << unique_config_elements.size() << std::endl;


  std::cout << "Prim config fg size: " << prim_config_fg.size() << std::endl;
  std::cout << "config scel fg size: " << config_scel_fg.size() << std::endl;


  /*
    int ii = 1;
    int kk = 4;
    SymOp Sik0 = prim_scel.prim_grid().sym_op(kk)*orbit.equivalence_map()[ii][0];

    int counter1 = 0;
    for(auto it : prim_config_fg) {
      SymRepIndexCompare compare;
      std::cout << counter1 << "\n";
      //it.sym_op().print(std::cout);
      //std::cout << "\n";
      counter1++;
      //std::cout << "Sik: " << std::endl;
      //Sik0.print(std::cout);
      //std::cout << "\n";
      SymOp itSik = it.sym_op()*Sik0;
      std::cout << "it.sym_op()*Sik: " << std::endl;
      itSik.print(std::cout);
      std::cout << "\n";
      std::cout << "Compared: " << compare(Sik0, it.sym_op()*Sik0) << std::endl;
    }




    /*
    std::vector<Kinetics::DiffusionTransformation> unique_config_difftrans;
    //Kinetics::sub_orbit_generators <Kinetics::DiffusionTransformation, Kinetics::PrimPeriodicDiffTransSymCompare, std::vector<Kinetics::DiffusionTransformation>::iterator> (orbit, config, std::back_inserter(unique_config_elements));
    Kinetics::sub_orbit_generators(orbit, config, std::back_inserter(unique_config_difftrans));

    std::cout << unique_config_difftrans.size() << std::endl;

    for(auto i : unique_config_difftrans) {
      std::cout << i << std::endl;
    }

    /*
    std::vector<Kinetics::DiffusionTransformation> unique_prim_elements;
    sub_orbit_generator(orbit, bg_config_prim, std::back_inserter(unique_prim_elements));


    for(i : orbit.size()) {
      for(k : scel.prim_grid.size()) {
        Sik = scel.prim_grid().sym_op(k)*orbit.equivalence_map()[i][0];

        auto lambda = [&](const PermuteIterator &it) {
          return Sik < it.sym_op*Sik;
        };

        if(std::none_of(config_scel_fg.begin(), config_scel_fg.end(), lambda)) {
          // use this Eik
        }
      }
    }

    /*

  */

  /*
  /// Find unique DiffusionTransformations
  //    PermuteIterator begin = bg_config_prim.supercell().permute_begin();
  //    PermuteIterator end = bg_config_prim.supercell().permute_end();
  //    Kinetics::DiffTransEnumEquivalents diff_trans_unique(diff_trans_prototype, begin, end, bg_config_prim);
  PermuteIterator begin = config.supercell().permute_begin();
  PermuteIterator end = config.supercell().permute_end();
  Kinetics::DiffTransEnumEquivalents diff_trans_unique(diff_trans_prototype, begin, end, config);

  //  std::cout << "Config factor groups" << "\n";
  //  for ( auto &g: diff_trans_unique.invariant_subgroup()){
  //  g.sym_op().print(std::cout);

  std::vector<Kinetics::DiffusionTransformation> subdifftrans;

  for(auto it = diff_trans_unique.begin(); it != diff_trans_unique.end(); ++it) {
    subdifftrans.push_back(*it);
  }
  /*std::cout << "Diff_Trans_Uniques " << subdifftrans.size() << " total" << "\n";
  for(auto it = subdifftrans.begin(); it != subdifftrans.end(); ++it) {
    std::cout << *it << "\n";
  }*/

  /// Check if number of unique diff trans is what we expect
  //Size of diff trans factor group should be equal to (size of scel fg)*(volume of scel)/(size of bg_config fg)
  /*std::cout << "size of Config factor group " << diff_trans_unique.invariant_subgroup().size() << std::endl;
  std::cout << "size of Supercell factor group " << bg_config_prim.supercell().factor_group().size() << std::endl;
  std::cout << "size of Supercell (#prims) " << bg_config_prim.supercell().volume() << std::endl;
  std::cout << "size of non-prim Supercell factor group " << config.supercell().factor_group().size() << std::endl;
  std::cout << "size of non-prim Supercell (#prims) " << config.supercell().volume() << std::endl;
  std::cout << "size of diff trans unique " << subdifftrans.size() << std::endl;
  */
  int config_fg_int = static_cast<int>(diff_trans_unique.invariant_subgroup().size());
  int scel_fg_int = static_cast<int>(bg_config_prim.supercell().factor_group().size());
  int scel_vol_int = static_cast<int>(bg_config_prim.supercell().volume());

  int subdifftrans_size_int = static_cast<int>(subdifftrans.size());

  int expected_subdifftrans_size_int = scel_fg_int * scel_vol_int ElementOutputIterator / config_fg_int;
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

  test::FCCTernaryProj proj2;
  proj2.check_init();
  proj2.check_composition();

  Logging logging2 = Logging::null();
  PrimClex primclex2(proj2.dir, logging2);

  std::vector<PrimPeriodicIntegralClusterOrbit> orbits2;
  make_prim_periodic_orbits(
    primclex2.prim(),
    bspecs,
    alloy_sites_filter,
    primclex2.crystallography_tol(),
    std::back_inserter(orbits2),
    primclex2.log());

  Eigen::Vector3d a2, b2, c2;
  std::tie(a2, b2, c2) = primclex2.prim().lattice().vectors();

  /// Make background_config
  Supercell fcc_scel {&primclex2, Lattice(2 * a2, 2 * b2, 2 * c2)};
  Configuration l12(fcc_scel);
  l12.init_occupation();
  l12.init_displacement();
  l12.init_deformation();
  l12.init_specie_id();
  l12.set_occupation({0, 0, 0, 1, 1, 0, 0, 0});

  /// Find prototype of m_diff_trans_orbit
  std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits2;
  Kinetics::make_prim_periodic_diff_trans_orbits(orbits2.begin() + 2, orbits2.begin() + 3, primclex2.crystallography_tol(), std::back_inserter(diff_trans_orbits2));
  Kinetics::DiffusionTransformation diff_trans_prototype2 = diff_trans_orbits2[4].prototype();

  //std::cout << "Prototype Diff Trans:" << "\n" << diff_trans_prototype2 << "\n";

  //For this hop and configuration
  // Expect 2 unique transformations A site to B site (16 total in scel) or A site to A site (32 total in scel)


  /// Make background_config
  Supercell fcc_scel2 {&primclex2, Lattice(4 * a2, 2 * b2, 2 * c2)};
  Configuration l12_ext(fcc_scel2);
  l12_ext.init_occupation();
  l12_ext.init_displacement();
  l12_ext.init_deformation();
  l12_ext.init_specie_id();
  l12_ext.set_occupation({0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0});

  ///For this hop and configuration
  ///Expect the A to B site along c axis (32 total) A to A site along a axis (32 total) and A to A along B (32 total)
}
BOOST_AUTO_TEST_SUITE_END()
