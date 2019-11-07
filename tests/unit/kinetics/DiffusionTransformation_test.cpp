#include "gtest/gtest.h"

/// What is being tested:
#include "casm/kinetics/DiffusionTransformation_impl.hh"
#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"

/// What is being used to test it:
#include "casm/app/AppIO.hh"
#include "casm/app/AppIO_impl.hh"
#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "ZrOProj.hh"
#include "casm/completer/Handlers.hh"
#include "casm/app/casm_functions.hh"
#include "casm/app/enum.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/casm_io/json_io/jsonFile.hh"
#include "casm/casm_io/container/stream_io.hh"

using namespace CASM;
using namespace test;

namespace {

  struct DiffTransEnumParserPtr {
    jsonParser input;
    std::unique_ptr<Kinetics::DiffTransEnumParser> parser;

    DiffTransEnumParserPtr(const PrimClex &primclex, const std::string &specs):
      input(jsonParser::parse(specs)),
      parser(notstd::make_unique<Kinetics::DiffTransEnumParser>(primclex, input, fs::path(), true)) {}

    Kinetics::DiffTransEnumParser *operator->() {
      return parser.get();
    }
  };
}

TEST(DiffusionTransformationTest, BasicsTest0) {

  typedef Kinetics::SpeciesLocation SpecieLocation;

  /// Make test project
  EXPECT_EQ(true, true);
  test::ZrOProj proj;
  proj.check_init();
  proj.check_composition();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);
  const Structure &prim = primclex.prim();

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();
  EXPECT_EQ(true, true);

  // Make background config
  Supercell _scel {&primclex, Lattice(1 * a, 1 * b, 1 * c)};
  Configuration _config(_scel);
  _config.set_occupation(std::vector<int>({0, 0, 1, 0}));

  Supercell background_scel {&primclex, Lattice(3 * a, 3 * b, 3 * c)};
  Configuration background_config = _config.
                                    fill_supercell(background_scel, primclex.prim().factor_group()).
                                    in_canonical_supercell();
  EXPECT_EQ(true, true);

  {
    // Construct
    Kinetics::DiffusionTransformation diff_trans(prim);
    EXPECT_EQ(true, true);
    EXPECT_EQ(diff_trans.occ_transform().size(), 0);

    UnitCellCoord uccoordA(prim, 3, 1, 1, 1);
    UnitCellCoord uccoordB(prim, 2, 1, 1, 2);
    Index iVa = 0;
    Index iO = 1;

    // Add transform (so that it's not sorted as constructed)
    diff_trans.occ_transform().emplace_back(uccoordB, iO, iVa);
    diff_trans.occ_transform().emplace_back(uccoordA, iVa, iO);
    EXPECT_EQ(true, true);
    EXPECT_EQ(diff_trans.is_valid_occ_transform(), true);

    // Add non-self consistent trajectory
    diff_trans.species_traj().emplace_back(SpecieLocation(uccoordA, iVa, 0), SpecieLocation(uccoordB, iVa, 0));
    diff_trans.species_traj().emplace_back(SpecieLocation(uccoordB, iVa, 0), SpecieLocation(uccoordA, iVa, 0));
    EXPECT_EQ(true, true);
    EXPECT_EQ(diff_trans.is_valid_species_traj(), true);
    EXPECT_EQ(diff_trans.species_types_map(), true);
    EXPECT_EQ(diff_trans.breaks_indivisible_mol(), false);
    EXPECT_EQ(diff_trans.is_subcluster_transformation(), false);
    EXPECT_EQ(diff_trans.is_self_consistent(), false);
    EXPECT_EQ(diff_trans.is_valid(), false);
    diff_trans.species_traj().clear();

    // Add valid trajectory
    diff_trans.species_traj().emplace_back(SpecieLocation(uccoordA, iVa, 0), SpecieLocation(uccoordB, iVa, 0));
    diff_trans.species_traj().emplace_back(SpecieLocation(uccoordB, iO, 0), SpecieLocation(uccoordA, iO, 0));
    EXPECT_EQ(true, true);
    EXPECT_EQ(diff_trans.is_valid_species_traj(), true);
    EXPECT_EQ(diff_trans.species_types_map(), true);
    EXPECT_EQ(diff_trans.breaks_indivisible_mol(), false);
    EXPECT_EQ(diff_trans.is_subcluster_transformation(), false);
    EXPECT_EQ(diff_trans.is_self_consistent(), true);
    EXPECT_EQ(diff_trans.is_valid(), true);

    // Check sorting
    EXPECT_EQ((diff_trans.occ_transform()[0].uccoord == uccoordA), false);
    EXPECT_EQ(diff_trans.is_sorted(), false);
    diff_trans.sort();
    EXPECT_EQ((diff_trans.occ_transform()[0].uccoord == uccoordA), true);
    EXPECT_EQ(true, true);
    EXPECT_EQ(diff_trans.is_sorted(), true);

    // Copy
    Kinetics::DiffusionTransformation other(diff_trans);

    // Compare
    EXPECT_EQ((other == diff_trans), true);

    // Check canonical form
    EXPECT_EQ(diff_trans.is_canonical(background_scel), false);

    {
      ScelIsCanonical<Kinetics::DiffusionTransformation> f(background_scel);
      EXPECT_EQ(f(diff_trans), false);
    }

    // Translate (to canonical form)
    other = diff_trans + UnitCell(1, 1, 1);

    // Compare
    EXPECT_EQ((diff_trans < other), true);

    // Check canonical form
    EXPECT_EQ(other.is_canonical(background_scel), true);
    EXPECT_EQ((diff_trans.canonical_form(background_scel) == other), true);

    {
      ScelIsCanonical<Kinetics::DiffusionTransformation> f(background_scel);
      EXPECT_EQ(f(other), true);
    }

    // Check invariant subgroup (haven't checked for accuracy)
    EXPECT_EQ(
      other.invariant_subgroup(
        prim.factor_group(),
        PrimPeriodicDiffTransSymCompare(primclex)).size(), 12);
    EXPECT_EQ(other.invariant_subgroup(background_scel).size(), 12);
  }

}

TEST(DiffusionTransformationTest, SpeedTest0) {

  test::ZrOProj proj;
  proj.check_init();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);

  fs::path bspecs_path = autotools::abs_srcdir() + "/tests/unit/kinetics/ZrO_bspecs_0.json";
  jsonParser bspecs {bspecs_path};

  std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
  make_prim_periodic_orbits(
    primclex.prim(),
    bspecs,
    alloy_sites_filter,
    primclex.crystallography_tol(),
    std::back_inserter(orbits),
    primclex.log());

  const auto &generating_grp = primclex.prim().factor_group();
  PrimPeriodicDiffTransSymCompare sym_compare {primclex.crystallography_tol()};

  OrbitGenerators<PrimPeriodicDiffTransOrbit> generators {generating_grp, sym_compare};

  for(auto orbit_it = orbits.begin(); orbit_it != orbits.end(); ++orbit_it) {
    Kinetics::DiffusionTransformationEnum e {orbit_it->prototype()};
    for(auto it = e.begin(); it != e.end(); ++it) {
      // generators.insert generates the sorted canonical form and inserts it
      generators.insert(*it);
    }
  }

  EXPECT_EQ(generators.elements.size(), 507);
}

TEST(DiffusionTransformationTest, SpeedTest1) {

  test::ZrOProj proj;
  proj.check_init();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);

  fs::path bspecs_path = autotools::abs_srcdir() + "/tests/unit/kinetics/ZrO_bspecs_0.json";
  jsonParser bspecs {bspecs_path};

  std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
  make_prim_periodic_orbits(
    primclex.prim(),
    bspecs,
    alloy_sites_filter,
    primclex.crystallography_tol(),
    std::back_inserter(orbits),
    primclex.log());

  // make orbits of DiffTrans
  std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
  Kinetics::make_prim_periodic_diff_trans_orbits(
    orbits.begin(),
    orbits.end(),
    primclex.crystallography_tol(),
    std::back_inserter(diff_trans_orbits),
    &primclex);
  EXPECT_EQ(true, true);

  EXPECT_EQ(diff_trans_orbits.size(), 507);
}


TEST(DiffusionTransformationTest, EnumTest0) {

  test::ZrOProj proj;
  proj.check_init();
  proj.check_composition();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);

  fs::path bspecs_path = autotools::abs_srcdir() + "/tests/unit/kinetics/ZrO_bspecs_0.json";
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
  EXPECT_EQ(trans.size(), 0);

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


  EXPECT_EQ(true, true);
  for(auto it = orbits.begin(); it != orbits.end(); ++it) {

    // print the IntegralCluster prototype used to generate DiffTrans
    // std::cout << "----------------\n";
    // std::cout << "IntegralCluster: \n" << it->prototype() << std::endl;
    // EXPECT_EQ(true, true);

    // check the number of valid DiffTrans
    Kinetics::DiffusionTransformationEnum e {it->prototype()};
    EXPECT_EQ(std::distance(e.begin(), e.end()), count_check[it->prototype().size()]);

    // make orbits of DiffTrans
    std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
    Kinetics::make_prim_periodic_diff_trans_orbits(
      it,
      it + 1,
      primclex.crystallography_tol(),
      std::back_inserter(diff_trans_orbits),
      &primclex);
    EXPECT_EQ(true, true);

    // check how many DiffTrans orbits there are for each IntegralCluster orbit
    orbit_count.push_back(diff_trans_orbits.size());
    if(expected_orbit_count_it != expected_orbit_count.end()) {
      EXPECT_EQ(diff_trans_orbits.size(), *expected_orbit_count_it);
      ++expected_orbit_count_it;
    }
    else {
      EXPECT_EQ(diff_trans_orbits.size(), 0);
    }

    // check the size of the DiffTrans orbits
    for(const auto &orb : diff_trans_orbits) {
      mult.push_back(orb.size());
      if(expected_mult_it != expected_mult.end()) {
        EXPECT_EQ(orb.size(), *expected_mult_it);
        ++expected_mult_it;
      }
      else {
        EXPECT_EQ(orb.size(), 0);
      }
    }
    jsonParser dtjson;
    dtjson.put_obj();
    if(diff_trans_orbits.size()) {
      to_json(diff_trans_orbits[0].prototype(), dtjson);
      Kinetics::DiffusionTransformation trans = jsonConstructor<Kinetics::DiffusionTransformation>::from_json(dtjson, primclex.prim());
      EXPECT_EQ(trans == diff_trans_orbits[0].prototype(), 1);
      EXPECT_EQ(diff_trans_orbits[0].prototype().is_valid(), 1);
      EXPECT_EQ(trans.is_valid(), 1);
    }

    //print DiffTrans prototypes
    //{
    //  PrototypePrinter<Kinetics::DiffusionTransformation> printer;
    //  print_clust(diff_trans_orbits.begin(), diff_trans_orbits.end(), std::cout, printer);
    //}
  }

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

TEST(DiffusionTransformationTest, DiffTransEnumParserTest) {
  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());

  {
    DiffTransEnumParserPtr parser(primclex, std::string(R"({
        "cspecs": {
          "orbit_branch_specs" : {
            "2" : {"max_length" : 5.0},
            "3" : {"max_length" : 3.0}
          }
        }
      })"));
    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->all_errors().size(), 0);
    EXPECT_EQ(parser->all_warnings().size(), 0);
    EXPECT_EQ(parser->cspecs().orbit_branch_specs().max_branch, 3);
    EXPECT_EQ(almost_equal(parser->cspecs().orbit_branch_specs().max_length(2), 5.), true);
    EXPECT_EQ(almost_equal(parser->cspecs().orbit_branch_specs().max_length(3), 3.), true);
    EXPECT_EQ(parser->cspecs().orbit_specs().custom_generators.elements.size(), 0);
    EXPECT_EQ(parser->dry_run(), false);
    EXPECT_EQ(parser->coord_type(), FRAC);
    EXPECT_EQ(parser->orbit_printer_opt().orbit_print_mode, ORBIT_PRINT_MODE::PROTO);
    EXPECT_EQ(parser->required_species().size(), 0);
    EXPECT_EQ(parser->excluded_species().size(), 0);

    primclex.log() << parser->report() << std::endl;
  }

  // missing cspecs
  {
    DiffTransEnumParserPtr parser(primclex, std::string(R"({
        "require": ["C"],
        "exclude": ["A", "B"],
        "dry_run": true
      })"));
    EXPECT_EQ(parser->valid(), false);
    EXPECT_EQ(parser->all_errors().size(), 1);
    EXPECT_EQ(parser->all_warnings().size(), 0);
    EXPECT_EQ(parser->cspecs().exists(), false);
    EXPECT_EQ(parser->dry_run(), true);
    EXPECT_EQ(parser->coord_type(), FRAC);
    EXPECT_EQ(parser->orbit_printer_opt().orbit_print_mode, ORBIT_PRINT_MODE::PROTO);
    EXPECT_EQ(parser->required_species().size(), 1);
    EXPECT_EQ(parser->excluded_species().size(), 2);

    primclex.log() << parser->report() << std::endl;
  }

  // missing cspecs data
  {
    DiffTransEnumParserPtr parser(primclex, std::string(R"({
        "cspecs" : {},
        "require": ["C"],
        "exclude": ["A", "B"],
        "dry_run": true
      })"));
    EXPECT_EQ(parser->valid(), false);
    EXPECT_EQ(parser->all_errors().size(), 1);
    EXPECT_EQ(parser->all_warnings().size(), 0);
    EXPECT_EQ(parser->cspecs().exists(), true);
    EXPECT_EQ(parser->cspecs().orbit_branch_specs().exists(), false);
    EXPECT_EQ(parser->cspecs().orbit_specs().exists(), false);
    EXPECT_EQ(parser->dry_run(), true);
    EXPECT_EQ(parser->coord_type(), FRAC);
    EXPECT_EQ(parser->orbit_printer_opt().orbit_print_mode, ORBIT_PRINT_MODE::PROTO);
    EXPECT_EQ(parser->required_species().size(), 1);
    EXPECT_EQ(parser->excluded_species().size(), 2);

    primclex.log() << parser->report() << std::endl;
  }

  // require and exclude
  {
    DiffTransEnumParserPtr parser(primclex, std::string(R"({
        "require": ["C"],
        "exclude": ["A", "B"],
        "cspecs": {
          "orbit_branch_specs" : {
            "2" : {"max_length" : 5.0},
            "3" : {"max_length" : 3.0}
          }
        },
        "dry_run": true
      })"));
    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->all_errors().size(), 0);
    EXPECT_EQ(parser->all_warnings().size(), 0);
    EXPECT_EQ(parser->required_species().size(), 1);
    EXPECT_EQ(parser->excluded_species().size(), 2);
    EXPECT_EQ(parser->dry_run(), true);

    primclex.log() << parser->report() << std::endl;
  }

  // "require" in wrong place
  {
    DiffTransEnumParserPtr parser(primclex, std::string(R"({
        "cspecs": {
          "require": ["Va"],
          "orbit_branch_specs" : {
            "2" : {"max_length" : 5.0},
            "3" : {"max_length" : 3.0}
          }
        },
        "dry_run": true
      })"));
    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->all_errors().size(), 0);
    EXPECT_EQ(parser->all_warnings().size(), 1);
    EXPECT_EQ(parser->dry_run(), true);

    primclex.log() << parser->report() << std::endl;
  }

  // "Va" not in prim
  {
    DiffTransEnumParserPtr parser(primclex, std::string(R"({
        "require": ["Va"],
        "cspecs": {
          "orbit_branch_specs" : {
            "2" : {"max_length" : 5.0},
            "3" : {"max_length" : 3.0}
          }
        },
        "dry_run": true
      })"));
    EXPECT_EQ(parser->valid(), false);
    EXPECT_EQ(parser->all_errors().size(), 1);
    EXPECT_EQ(parser->all_warnings().size(), 0);
    EXPECT_EQ(parser->required_species().size(), 1);
    EXPECT_EQ(parser->dry_run(), true);

    primclex.log() << parser->report() << std::endl;
  }

  // "unrecognized" option
  {
    DiffTransEnumParserPtr parser(primclex, std::string(R"({
        "cspecs": {
          "orbit_branch_specs" : {
            "2" : {"max_length" : 5.0},
            "3" : {"max_length" : 3.0}
          }
        },
        "dry_run": true,
        "unrecognized": true
      })"));
    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->all_errors().size(), 0);
    EXPECT_EQ(parser->all_warnings().size(), 1);
    EXPECT_EQ(parser->dry_run(), true);

    primclex.log() << parser->report() << std::endl;
  }

  // bad "coordinate_mode" value
  {
    DiffTransEnumParserPtr parser(primclex, std::string(R"({
        "cspecs": {
          "orbit_branch_specs" : {
            "2" : {"max_length" : 5.0},
            "3" : {"max_length" : 3.0}
          }
        },
        "dry_run": true,
        "coordinate_mode": true
      })"));
    EXPECT_EQ(parser->valid(), false);
    EXPECT_EQ(parser->all_errors().size(), 1);
    EXPECT_EQ(parser->all_warnings().size(), 0);
    EXPECT_EQ(parser->dry_run(), true);

    primclex.log() << parser->report() << std::endl;
  }
}

TEST(DiffusionTransformationTest, EnumTest1) {

  test::ZrOProj proj;
  proj.check_init();
  proj.check_composition();

  auto &log = default_log();
  log.set_verbosity(Log::verbose);
  PrimClex primclex(proj.dir, log);

  log.custom("Prim");
  jsonParser tjson;
  write_prim(primclex.prim(), tjson, FRAC);
  log << tjson << std::endl << std::endl;
  write_prim(primclex.prim(), tjson, CART);
  log << tjson << std::endl;

  SymInfoOptions opt;

  opt.coord_type = FRAC;
  log.custom("Prim factor group (FRAC)");
  description(log, primclex.prim().factor_group(), primclex.prim().lattice(), opt);
  log << std::endl;

  opt.coord_type = CART;
  log.custom("Prim factor group (CART)");
  description(log, primclex.prim().factor_group(), primclex.prim().lattice(), opt);
  log << std::endl;

  opt.prec = 5;
  opt.coord_type = FRAC;
  log.custom("Prim factor group (brief FRAC)");
  brief_description(log, primclex.prim().factor_group(), primclex.prim().lattice(), opt);
  log << std::endl;

  opt.coord_type = CART;
  log.custom("Prim factor group (brief CART)");
  brief_description(log, primclex.prim().factor_group(), primclex.prim().lattice(), opt);
  log << std::endl;

  log.custom(Kinetics::DiffusionTransformationEnum::enumerator_name + " interface help");
  log << Kinetics::DiffusionTransformationEnum::interface_help() << std::endl;

  log.custom("OrbitPrinterOptionsParser standard help");
  log << OrbitPrinterOptionsParser::standard_help() << std::endl;

  log.custom("SymInfoOptionsParser interface help");
  log << SymInfoOptionsParser::standard_help() << std::endl;


  jsonFile diff_trans_json {autotools::abs_srcdir() + "/tests/unit/kinetics/ZrO_diff_trans_0.json"};
  diff_trans_json["dry_run"] = true;
  diff_trans_json["coordinate_mode"] = std::string("CART");
  diff_trans_json["orbit_printer_opt"]["print_invariant_grp"] = true;
  std::cout << "\ndiff_trans_json:\n" << diff_trans_json << std::endl;
  Completer::EnumOption enum_opt;
  enum_opt.desc();
  int success = Kinetics::DiffusionTransformationEnum::run(primclex, diff_trans_json, enum_opt, nullptr);
  EXPECT_EQ(primclex.generic_db<Kinetics::PrimPeriodicDiffTransOrbit>().size(), 123);
  EXPECT_EQ(success, 0);
}

TEST(DiffusionTransformationTest, EnumTest2) {

  test::FCCTernaryProj proj;
  proj.check_init();
  proj.check_composition();

  auto &log = null_log();
  log.set_verbosity(Log::verbose);
  PrimClex primclex(proj.dir, log);

  log.custom("Prim");
  jsonParser tjson;
  write_prim(primclex.prim(), tjson, FRAC);
  log << tjson << std::endl << std::endl;
  write_prim(primclex.prim(), tjson, CART);
  log << tjson << std::endl;

  SymInfoOptions opt;

  opt.coord_type = FRAC;
  log.custom("Prim factor group (FRAC)");
  description(log, primclex.prim().factor_group(), primclex.prim().lattice(), opt);
  log << std::endl;

  opt.coord_type = CART;
  log.custom("Prim factor group (CART)");
  description(log, primclex.prim().factor_group(), primclex.prim().lattice(), opt);
  log << std::endl;

  opt.prec = 5;
  opt.coord_type = FRAC;
  log.custom("Prim factor group (brief FRAC)");
  brief_description(log, primclex.prim().factor_group(), primclex.prim().lattice(), opt);
  log << std::endl;

  opt.coord_type = CART;
  log.custom("Prim factor group (brief CART)");
  brief_description(log, primclex.prim().factor_group(), primclex.prim().lattice(), opt);
  log << std::endl;


  jsonFile diff_trans_json {autotools::abs_srcdir() + "/tests/unit/kinetics/FCCTernary_diff_trans_err_0.json"};
  diff_trans_json["coordinate_mode"] = std::string("FRAC");
  diff_trans_json["orbit_printer_opt"]["orbit_print_mode"] = std::string("FULL");
  diff_trans_json["orbit_printer_opt"]["print_invariant_grp"] = true;
  diff_trans_json["orbit_printer_opt"]["print_equivalence_map"] = true;
  diff_trans_json["orbit_printer_opt"]["sym_info_opt"]["prec"] = 3;
  diff_trans_json["orbit_printer_opt"]["sym_info_opt"]["print_matrix_tau"] = true;
  Completer::EnumOption enum_opt;
  enum_opt.desc();
  int success = Kinetics::DiffusionTransformationEnum::run(primclex, diff_trans_json, enum_opt, nullptr);
  EXPECT_EQ(primclex.generic_db<Kinetics::PrimPeriodicDiffTransOrbit>().size(), 28);
  EXPECT_EQ(success, 1);
}
