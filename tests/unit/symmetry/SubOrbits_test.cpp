#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/symmetry/SubOrbits.hh"
#include "casm/symmetry/SubOrbits_impl.hh"

#include "casm/symmetry/ScelSubOrbits.hh"
#include "casm/symmetry/ScelSubOrbits_impl.hh"

#include "casm/symmetry/ConfigSubOrbits.hh"
#include "casm/symmetry/ConfigSubOrbits_impl.hh"

/// What is being used to test it:
#include "Common.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"

#include "casm/casm_io/VaspIO.hh"


using namespace CASM;

BOOST_AUTO_TEST_SUITE(SubOrbitsTest)

/// Test individual parts of MakeConfigSubOrbitGenerators
BOOST_AUTO_TEST_CASE(ZrOProj) {

  /// Make test project
  BOOST_CHECK_EQUAL(true, true);
  test::ZrOProj proj;
  proj.check_init();
  proj.check_composition();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);
  const Structure &prim = primclex.prim();
  const Lattice &lat = prim.lattice();
  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();
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
  BOOST_CHECK_EQUAL(orbits.size(), 74); // value checked with casm0.2x

  //print_clust(orbits.begin(), orbits.end(), std::cout, PrototypePrinter<IntegralCluster>());

  // Make PrimPeriodicDiffTransOrbit
  std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
  Kinetics::make_prim_periodic_diff_trans_orbits(
    orbits.begin() + 2,
    orbits.begin() + 4,
    primclex.crystallography_tol(),
    std::back_inserter(diff_trans_orbits),
    &primclex);
  BOOST_CHECK_EQUAL(true, true);
  BOOST_CHECK_EQUAL(diff_trans_orbits.size(), 4);
  // O-O pair within layer and O-O pair between layers
  // O-O and O-Va for each case -> 4 total diff_trans_orbits

  /*
  std::cout << "diff_trans_orbits: \n" << std::endl;
  print_clust(
    diff_trans_orbits.begin(),
    diff_trans_orbits.end(),
    std::cout,
    PrototypePrinter<Kinetics::DiffusionTransformation>());
  */

  Printer<Kinetics::DiffusionTransformation> diff_trans_printer;
  typedef std::vector<Kinetics::DiffusionTransformation> DiffTransVec;
  typedef ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> ScelDiffTransSymCompare;

  // Test 0 Step-by-step of internals of MakeConfigSubOrbitGenerators
  {
    /*
    Supercell scel {&primclex, Lattice(2 * a, 1 * b, 1 * c)};
    Configuration config(scel);
    config.set_occupation({0, 0, 0, 0, 1, 1, 0, 0});
    */

    Supercell scel {&primclex, Lattice(3 * a, 2 * b, 1 * c)};
    Configuration config(scel);
    config.set_occupation({
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      1, 1, 1, 1, 1, 1,
      0, 0, 0, 0, 0, 0
    });

    {
      std::vector<IntegralCluster> generators;
      BOOST_CHECK_EQUAL(true, true);
      make_suborbit_generators(orbits[0], config, std::back_inserter(generators));
      //std::cout << "orbit 0: \n" << orbits[0].prototype() << std::endl;
      BOOST_CHECK_EQUAL(generators.size(), 1);
    }

    {
      std::vector<IntegralCluster> generators;
      BOOST_CHECK_EQUAL(true, true);
      make_suborbit_generators(orbits[1], config, std::back_inserter(generators));
      //std::cout << "orbit 1: \n" << orbits[1].prototype() << std::endl;
      BOOST_CHECK_EQUAL(generators.size(), 2);
    }

    {
      std::vector<IntegralCluster> generators;
      BOOST_CHECK_EQUAL(true, true);
      make_suborbit_generators(orbits[2], config, std::back_inserter(generators));
      //std::cout << "orbit 2: \n" << orbits[2].prototype() << std::endl;
      BOOST_CHECK_EQUAL(generators.size(), 1);
    }

    ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> scel_sym_compare(
      config.supercell().prim_grid(),
      config.crystallography_tol());
    // Configuration with every other layer of O filled


    Configuration prim_config = config.primitive().in_canonical_supercell();
    std::vector<PermuteIterator> prim_config_fg = prim_config.factor_group();
    ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> prim_sym_compare(
      prim_config.supercell().prim_grid(),
      prim_config.crystallography_tol());

    BOOST_CHECK_EQUAL(true, true);

    //std::cout << "config: \n" << config << std::endl;
    //std::cout << "prim_config: \n" << prim_config << std::endl;

    auto print_res = [&](std::string name, const DiffTransVec & diff_trans_vec, ScelDiffTransSymCompare sym_compare) {
      std::cout << name << ":" << std::endl;
      Index index = 0;
      for(const auto &el : diff_trans_vec) {
        std::cout << "index: " << index << std::endl;
        diff_trans_printer.print(sym_compare.prepare(el), std::cout);
        ++index;
      }
      std::cout << std::endl;
    };

    /// prim -> prim_config.supercell() symmetry breaking
    DiffTransVec prim_scel_suborbit_generators;
    {
      MakeSubOrbitGenerators gen(
        prim_config.prim().factor_group(),
        prim_config.supercell().factor_group());
      for(const auto &orbit : diff_trans_orbits) {
        gen(orbit, std::back_inserter(prim_scel_suborbit_generators));
      }
      BOOST_CHECK_EQUAL(true, true);
      BOOST_CHECK_EQUAL(prim_scel_suborbit_generators.size(), 4);
      // prim and prim_config.supercell() are same size and shape
      // O-Va hop in c direction, O-O hop in c direction
      // O-Va hop in a direction, O-O hop in a direction; a = b
      // total of 4 distinct hops

      //print_res("prim_scel_suborbit_generators", prim_scel_suborbit_generators, prim_sym_compare);
    }

    /// prim_config.supercell() -> prim_config symmetry breaking
    std::vector<Kinetics::DiffusionTransformation> prim_config_suborbit_generators;
    {
      for(const auto &el : prim_scel_suborbit_generators) {
        make_suborbit_generators(
          el,
          prim_config.supercell(),
          prim_config_fg.begin(),
          prim_config_fg.end(),
          std::back_inserter(prim_config_suborbit_generators));
      }
      BOOST_CHECK_EQUAL(true, true);
      BOOST_CHECK_EQUAL(prim_config_suborbit_generators.size(), 6);
      // O-Va hop in c direction, O-O hop in c direction
      // In Va layer: O-Va hop in a direction, O-O hop in a direction; a = b
      // In O layer: O-Va hop in a direction, O-O hop in a direction; a = b
      // total of 6 distinct hops

      //print_res("prim_config_suborbit_generators", prim_config_suborbit_generators, prim_sym_compare);
    }

    /// prim_config -> config symmetry breaking
    std::vector<Kinetics::DiffusionTransformation> config_suborbit_generators;
    {

      //std::cout << "prim_config.supercell().factor_group().size(): "
      //  << prim_config.supercell().factor_group().size() << std::endl;
      //std::cout << "config.supercell().factor_group().size(): "
      //  << config.supercell().factor_group().size() << std::endl;

      std::vector<PermuteIterator> config_subgroup = make_invariant_subgroup(
                                                       config.supercell(),
                                                       prim_config.supercell(),
                                                       prim_config_fg.begin(),
                                                       prim_config_fg.end());

      for(const auto &el : prim_config_suborbit_generators) {
        make_suborbit_generators(
          el,
          prim_config.supercell(),
          prim_config_fg.begin(), prim_config_fg.end(),
          config_subgroup.begin(), config_subgroup.end(),
          std::back_inserter(config_suborbit_generators));
      }
      BOOST_CHECK_EQUAL(true, true);
      BOOST_CHECK_EQUAL(config_suborbit_generators.size(), 14);
      // O-Va hop in c direction, O-O hop in c direction
      // In Va layer: O-Va hop in a/b/a-b directions, O-O hop in a/b/a-b directions; a -=! b =! a-b
      // In O layer: O-Va hop in a/b/a-b directions, O-O hop in a/b/a-b directions; a -=! b =! a-b
      // total of 14 distinct hops

      //print_res("config_suborbit_generators 0", config_suborbit_generators, scel_sym_compare);
    }

    // test all combined in MakeConfigSubOrbitGenerators
    {
      std::vector<Kinetics::DiffusionTransformation> config_suborbit_generators;
      MakeConfigSubOrbitGenerators{config}(
        diff_trans_orbits.begin(),
        diff_trans_orbits.end(),
        std::back_inserter(config_suborbit_generators));
      BOOST_CHECK_EQUAL(true, true);
      BOOST_CHECK_EQUAL(config_suborbit_generators.size(), 14);
      // checked against prim_config -> config test

      //print_res("config_suborbit_generators 1", config_suborbit_generators, scel_sym_compare);
    }

    // test all combined in MakeConfigSubOrbitGenerators
    {
      std::vector<Kinetics::DiffusionTransformation> config_suborbit_generators;
      make_suborbit_generators(
        diff_trans_orbits.begin(),
        diff_trans_orbits.end(),
        config,
        std::back_inserter(config_suborbit_generators));
      BOOST_CHECK_EQUAL(true, true);
      BOOST_CHECK_EQUAL(config_suborbit_generators.size(), 14);
      // checked against prim_config -> config test

      //print_res("config_suborbit_generators 2", config_suborbit_generators, scel_sym_compare);
    }

  }
}

/// FCC structure with pure species A background
/// and nearest-neighbor hop species A, B, and C
BOOST_AUTO_TEST_CASE(FCCTernaryProj) {

  /// Make test project
  BOOST_CHECK_EQUAL(true, true);
  test::FCCTernaryProj proj;
  proj.check_init();
  proj.check_composition();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);
  const Structure &prim = primclex.prim();
  const Lattice &lat = prim.lattice();
  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();
  BOOST_CHECK_EQUAL(true, true);

  // Make PrimPeriodicIntegralClusterOrbit
  fs::path bspecs_path = "tests/unit/kinetics/bspecs_1.json";
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
  BOOST_CHECK_EQUAL(orbits.size(), 31); // value checked with casm0.2x

  //print_clust(orbits.begin(), orbits.end(), std::cout, PrototypePrinter<IntegralCluster>());

  // Make PrimPeriodicDiffTransOrbit
  std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
  Kinetics::make_prim_periodic_diff_trans_orbits(
    orbits.begin() + 2,
    orbits.begin() + 3,
    primclex.crystallography_tol(),
    std::back_inserter(diff_trans_orbits),
    &primclex);
  BOOST_CHECK_EQUAL(true, true);
  BOOST_CHECK_EQUAL(diff_trans_orbits.size(), 6);
  // nearest neighbor A-A, A-B, A-C, B-B, B-C, C-C pairs; total of 6

  /*
  std::cout << "diff_trans_orbits: \n" << std::endl;
  print_clust(
    diff_trans_orbits.begin(),
    diff_trans_orbits.end(),
    std::cout,
    PrototypePrinter<Kinetics::DiffusionTransformation>());
  /**/

  Printer<Kinetics::DiffusionTransformation> diff_trans_printer;
  typedef std::vector<Kinetics::DiffusionTransformation> DiffTransVec;
  typedef ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> ScelDiffTransSymCompare;

  // Test 0 Step-by-step of internals of MakeConfigSubOrbitGenerators
  {
    /*
    Supercell scel {&primclex, Lattice(2 * a, 1 * b, 1 * c)};
    Configuration config(scel);
    config.set_occupation({0, 0, 0, 0, 1, 1, 0, 0});
    */

    Supercell scel {&primclex, Lattice(3 * a, 2 * b, 2 * c)};
    Configuration config(scel);
    config.set_occupation({
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0
    });
    ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> scel_sym_compare(
      config.supercell().prim_grid(),
      config.crystallography_tol());

    // Print out a POSCAR file
    // fs::ofstream file;
    //   fs::path POSCARpath = "POSCAR";
    //   file.open(POSCARpath);
    //   VaspIO::PrintPOSCAR p(config);
    //   p.sort();
    //   p.print(file);
    //   file.close();

    Configuration prim_config = config.primitive().in_canonical_supercell();
    std::vector<PermuteIterator> prim_config_fg = prim_config.factor_group();
    ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> prim_sym_compare(
      prim_config.supercell().prim_grid(),
      prim_config.crystallography_tol());

    BOOST_CHECK_EQUAL(true, true);

    //std::cout << "config: \n" << config << std::endl;
    //std::cout << "prim_config: \n" << prim_config << std::endl;

    auto print_res = [&](std::string name, const DiffTransVec & diff_trans_vec, ScelDiffTransSymCompare sym_compare) {
      std::cout << name << ":" << std::endl;
      Index index = 0;
      for(const auto &el : diff_trans_vec) {
        std::cout << "index: " << index << std::endl;
        diff_trans_printer.print(sym_compare.prepare(el), std::cout);
        ++index;
      }
      std::cout << std::endl;
    };

    /// prim -> prim_config.supercell() symmetry breaking
    DiffTransVec prim_scel_suborbit_generators;
    {
      MakeSubOrbitGenerators gen(
        prim_config.prim().factor_group(),
        prim_config.supercell().factor_group());
      for(const auto &orbit : diff_trans_orbits) {
        gen(orbit, std::back_inserter(prim_scel_suborbit_generators));
      }
      BOOST_CHECK_EQUAL(true, true);
      BOOST_CHECK_EQUAL(prim_scel_suborbit_generators.size(), 6);
      // prim and prim_config.supercell() are equivalent
      // nearest neighbor A-A, A-B, A-C, B-B, B-C, C-C pairs; total of 6

      //print_res("prim_scel_suborbit_generators", prim_scel_suborbit_generators, prim_sym_compare);
    }

    /// prim_config.supercell() -> prim_config symmetry breaking
    std::vector<Kinetics::DiffusionTransformation> prim_config_suborbit_generators;
    {
      for(const auto &el : prim_scel_suborbit_generators) {
        make_suborbit_generators(
          el,
          prim_config.supercell(),
          prim_config_fg.begin(),
          prim_config_fg.end(),
          std::back_inserter(prim_config_suborbit_generators));
      }
      BOOST_CHECK_EQUAL(true, true);
      BOOST_CHECK_EQUAL(prim_config_suborbit_generators.size(), 6);
      // prim_config.supercell() and prim_config are equivalent
      // nearest neighbor A-A, A-B, A-C, B-B, B-C, C-C pairs; total of 6

      //print_res("prim_config_suborbit_generators", prim_config_suborbit_generators, prim_sym_compare);
    }

    /// prim_config -> config symmetry breaking
    std::vector<Kinetics::DiffusionTransformation> config_suborbit_generators;
    {

      // std::cout << "prim_config.supercell().factor_group().size(): "
      //  << prim_config.supercell().factor_group().size() << std::endl;
      // std::cout << "config.supercell().factor_group().size(): "
      //  << config.supercell().factor_group().size() << std::endl;

      std::vector<PermuteIterator> config_subgroup = make_invariant_subgroup(
                                                       config.supercell(),
                                                       prim_config.supercell(),
                                                       prim_config_fg.begin(),
                                                       prim_config_fg.end());

      for(const auto &el : prim_config_suborbit_generators) {
        make_suborbit_generators(
          el,
          prim_config.supercell(),
          prim_config_fg.begin(), prim_config_fg.end(),
          config_subgroup.begin(), config_subgroup.end(),
          std::back_inserter(config_suborbit_generators));
      }
      BOOST_CHECK_EQUAL(true, true);
      BOOST_CHECK_EQUAL(config_suborbit_generators.size(), 24);
      // b-c, a, a-b, and c directions; b = c != a

      //print_res("config_suborbit_generators 0", config_suborbit_generators, scel_sym_compare);
    }

    // test all combined in MakeConfigSubOrbitGenerators
    {
      std::vector<Kinetics::DiffusionTransformation> config_suborbit_generators;
      MakeConfigSubOrbitGenerators{config}(
        diff_trans_orbits.begin(),
        diff_trans_orbits.end(),
        std::back_inserter(config_suborbit_generators));
      BOOST_CHECK_EQUAL(true, true);
      BOOST_CHECK_EQUAL(config_suborbit_generators.size(), 24);

      //print_res("config_suborbit_generators 1", config_suborbit_generators, scel_sym_compare);
      // checked against prim_config -> config test
    }

    // test all combined in MakeConfigSubOrbitGenerators
    {
      std::vector<Kinetics::DiffusionTransformation> config_suborbit_generators;
      make_suborbit_generators(
        diff_trans_orbits.begin(),
        diff_trans_orbits.end(),
        config,
        std::back_inserter(config_suborbit_generators));
      BOOST_CHECK_EQUAL(true, true);
      BOOST_CHECK_EQUAL(config_suborbit_generators.size(), 24);

      //print_res("config_suborbit_generators 2", config_suborbit_generators, scel_sym_compare);
      // checked against prim_config -> config test
    }

  }
}

/// FCC structure with L12 ordering
/// and nearest-neighbor ternary hop of species A
BOOST_AUTO_TEST_CASE(L12Proj) {

  /// Make test project
  BOOST_CHECK_EQUAL(true, true);
  test::FCCTernaryProj proj;
  proj.check_init();
  proj.check_composition();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);
  const Structure &prim = primclex.prim();
  const Lattice &lat = prim.lattice();
  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();
  BOOST_CHECK_EQUAL(true, true);

  // Make PrimPeriodicIntegralClusterOrbit
  fs::path bspecs_path = "tests/unit/kinetics/bspecs_1.json";
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
  BOOST_CHECK_EQUAL(orbits.size(), 31); // same as FCC

  // Make FCC standard unit cell
  Eigen::Matrix3i superlattice_matrix;
  superlattice_matrix <<
                      -1, 1, 1,
                      1, -1, 1,
                      1, 1, -1;
  Supercell scel {&primclex, superlattice_matrix};

  //print_clust(orbits.begin(), orbits.end(), std::cout, PrototypePrinter<IntegralCluster>());

  // Make PrimPeriodicDiffTransOrbit
  std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
  Kinetics::make_prim_periodic_diff_trans_orbits(
    orbits.begin() + 8,
    orbits.begin() + 9,
    primclex.crystallography_tol(),
    std::back_inserter(diff_trans_orbits),
    &primclex);
  BOOST_CHECK_EQUAL(true, true);
  BOOST_CHECK_EQUAL(diff_trans_orbits.size(), 10);
  // Nearest neighbor triplet
  // A-A-A, A-A-B, A-A-C, A-B-B, A-B-C, A-C-C
  // B-B-B, B-B-C, B-C-C, C-C-C
  // total of 10 unique orbits

  /*
  std::cout << "diff_trans_orbits: \n" << std::endl;
  print_clust(
    diff_trans_orbits.begin(),
    diff_trans_orbits.end(),
    std::cout,
    PrototypePrinter<Kinetics::DiffusionTransformation>());
  /**/

  Printer<Kinetics::DiffusionTransformation> diff_trans_printer;
  typedef std::vector<Kinetics::DiffusionTransformation> DiffTransVec;
  typedef ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> ScelDiffTransSymCompare;

  // Test 0 Step-by-step of internals of MakeConfigSubOrbitGenerators
  {

    // Decorating FCC to have L12
    //Supercell scel {&primclex, Lattice(1 * a, 1 * b, 1 * c)};
    Configuration config(scel);
    config.set_occupation({
      0, 1, 1, 1
    });

    // supercell of standard FCC unit cell
    Eigen::Matrix3i superlattice_matrix_2;
    superlattice_matrix_2 <<
                          -1, 1, 2,
                          1, -1, 2,
                          1, 1, -2;
    Supercell bg_scel {&primclex, superlattice_matrix_2};
    Configuration bg_config = config.fill_supercell(bg_scel, primclex.prim().factor_group());

    ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> scel_sym_compare(
      config.supercell().prim_grid(),
      config.crystallography_tol());


    Configuration prim_config = bg_config.primitive().in_canonical_supercell();
    std::vector<PermuteIterator> prim_config_fg = prim_config.factor_group();
    ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> prim_sym_compare(
      prim_config.supercell().prim_grid(),
      prim_config.crystallography_tol());

    BOOST_CHECK_EQUAL(true, true);

    /// Print out a POSCAR file
    // fs::ofstream file;
    //   fs::path POSCARpath = "POSCAR";
    //   file.open(POSCARpath);
    //   VaspIO::PrintPOSCAR p(bg_config);
    //   p.sort();
    //   p.print(file);
    //   file.close();

    // std::cout << "config: \n" << config << std::endl;
    // std::cout << "prim_config: \n" << prim_config << std::endl;

    auto print_res = [&](std::string name, const DiffTransVec & diff_trans_vec, ScelDiffTransSymCompare sym_compare) {
      std::cout << name << ":" << std::endl;
      Index index = 0;
      for(const auto &el : diff_trans_vec) {
        std::cout << "index: " << index << std::endl;
        diff_trans_printer.print(sym_compare.prepare(el), std::cout);
        ++index;
      }
      std::cout << std::endl;
    };

    /// prim -> prim_config.supercell() symmetry breaking
    DiffTransVec prim_scel_suborbit_generators;
    {
      MakeSubOrbitGenerators gen(
        prim_config.prim().factor_group(),
        prim_config.supercell().factor_group());
      for(const auto &orbit : diff_trans_orbits) {
        gen(orbit, std::back_inserter(prim_scel_suborbit_generators));
      }
      BOOST_CHECK_EQUAL(true, true);
      BOOST_CHECK_EQUAL(prim_scel_suborbit_generators.size(), 10);
      // prim is equivalent to prim_config.supercell()

      //print_res("prim_scel_suborbit_generators", prim_scel_suborbit_generators, prim_sym_compare);
    }

    /// prim_config.supercell() -> prim_config symmetry breaking
    std::vector<Kinetics::DiffusionTransformation> prim_config_suborbit_generators;
    // {
    //   for(const auto &el : prim_scel_suborbit_generators) {
    //     make_suborbit_generators(
    //       el,
    //       prim_config.supercell(),
    //       prim_config_fg.begin(),
    //       prim_config_fg.end(),
    //       std::back_inserter(prim_config_suborbit_generators));
    //   }
    // Use B-A-B cluster
    {
      make_suborbit_generators(
        prim_scel_suborbit_generators[6],
        prim_config.supercell(),
        prim_config_fg.begin(),
        prim_config_fg.end(),
        std::back_inserter(prim_config_suborbit_generators));

      BOOST_CHECK_EQUAL(true, true);
      BOOST_CHECK_EQUAL(prim_config_suborbit_generators.size(), 3);
      // Rotation on cluster with all B sites
      // Rotation on 1A/2B cluster sites with A going from A site to B site
      // Rotation on 1A/2B cluster sites with A going from B site to B site

      //print_res("prim_config_suborbit_generators", prim_config_suborbit_generators, prim_sym_compare);
    }

    /// prim_config -> config symmetry breaking
    std::vector<Kinetics::DiffusionTransformation> config_suborbit_generators;
    {

      // std::cout << "prim_config.supercell().factor_group().size(): "
      //  << prim_config.supercell().factor_group().size() << std::endl;
      // std::cout << "bg_config.supercell().factor_group().size(): "
      //  << bg_config.supercell().factor_group().size() << std::endl;

      std::vector<PermuteIterator> config_subgroup = make_invariant_subgroup(
                                                       bg_config.supercell(),
                                                       prim_config.supercell(),
                                                       prim_config_fg.begin(),
                                                       prim_config_fg.end());

      for(const auto &el : prim_config_suborbit_generators) {
        make_suborbit_generators(
          el,
          prim_config.supercell(),
          prim_config_fg.begin(), prim_config_fg.end(),
          config_subgroup.begin(), config_subgroup.end(),
          std::back_inserter(config_suborbit_generators));
      }
      BOOST_CHECK_EQUAL(true, true);
      BOOST_CHECK_EQUAL(config_suborbit_generators.size(), 7);
      // LB = B site on long side; SB = B site on short site; A = A site
      // Cluster of LB-SB-A with A species moving from SB to LB site
      // Cluster of LB-LB-A with A species moving from LB to LB site
      // Cluster of SB-LB-LB with A species moving from LB to SB site
      // Cluster of SB-LB-LB with A species moving from LB to LB site
      // Cluster of LB-LB-A with A species moving from LB to A site
      // Cluster of LB-SB-A with A species moving from A to LB site
      // Cluster of LB-SB-A with A species moving from A to SB site
      // total of 7 unique diff_trans

      //print_res("config_suborbit_generators 0", config_suborbit_generators, scel_sym_compare);
    }

    // test all combined in MakeConfigSubOrbitGenerators
    {
      std::vector<Kinetics::DiffusionTransformation> config_suborbit_generators;
      MakeConfigSubOrbitGenerators{bg_config}(
        diff_trans_orbits[6],
        std::back_inserter(config_suborbit_generators));
      BOOST_CHECK_EQUAL(true, true);
      BOOST_CHECK_EQUAL(config_suborbit_generators.size(), 7);
      // checked against prim_config -> config test

      //print_res("config_suborbit_generators 1", config_suborbit_generators, scel_sym_compare);
    }

    // test all combined in MakeConfigSubOrbitGenerators
    {
      std::vector<Kinetics::DiffusionTransformation> config_suborbit_generators;
      make_suborbit_generators(
        diff_trans_orbits[6],
        bg_config,
        std::back_inserter(config_suborbit_generators));
      BOOST_CHECK_EQUAL(true, true);
      BOOST_CHECK_EQUAL(config_suborbit_generators.size(), 7);
      // checked against prim_config -> config test

      //print_res("config_suborbit_generators 2", config_suborbit_generators, scel_sym_compare);
    }

  }
}

BOOST_AUTO_TEST_SUITE_END()
