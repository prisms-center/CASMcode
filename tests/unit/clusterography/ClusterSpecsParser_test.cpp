#include "gtest/gtest.h"
#include <memory>

/// What is being tested:
#include "casm/clex/Supercell.hh"
#include "casm/clusterography/ClusterSpecsParser_impl.hh"

/// What is being used to test it:
#include "FCCTernaryProj.hh"
#include "TestConfiguration.hh"
#include "TestOrbits.hh"
#include "TestSupercell.hh"
#include "ZrOProj.hh"
#include "casm/casm_io/json/jsonFile.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/completer/Handlers.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"
#include "casm/kinetics/DiffusionTransformation_impl.hh"
#include "casm/misc/CRTPBase.hh"
#include "casm/symmetry/ElementSymApply.hh"

using namespace CASM;

namespace {

  struct PrimParser {
    jsonParser input;
    std::unique_ptr<PrimPeriodicClustersByMaxLength> parser;

    PrimParser(const PrimClex &primclex, const std::string &specs):
      input(jsonParser::parse(specs)),
      parser(notstd::make_unique<PrimPeriodicClustersByMaxLength>(
               primclex,
               primclex.prim().factor_group(),
               PrimPeriodicSymCompare<IntegralCluster>(primclex.shared_prim(), primclex.crystallography_tol()),
               input,
               fs::path(),
               true)) {}
    PrimPeriodicClustersByMaxLength *operator->() {
      return parser.get();
    }
  };

  struct LocalParser {
    jsonParser input;
    std::unique_ptr<LocalClustersByMaxLength<Kinetics::DiffusionTransformation>> parser;

    LocalParser(const Supercell &scel, const std::string &specs):
      input(jsonParser::parse(specs)),
      parser(notstd::make_unique<LocalClustersByMaxLength<Kinetics::DiffusionTransformation>>(
               scel,
               input,
               fs::path(),
               true)) {}
    LocalClustersByMaxLength<Kinetics::DiffusionTransformation> *operator->() {
      return parser.get();
    }
  };

  struct TestOrbits0 : test::TestPrimPeriodicDiffusionTransformationOrbits {
    TestOrbits0(const PrimClex &primclex) :
      test::TestPrimPeriodicDiffusionTransformationOrbits(
        primclex,
        jsonFile(autotools::abs_srcdir() + "/tests/unit/kinetics/FCCTernary_bspecs_0.json"),
        2, 4) {
      // counts only check for consistency, not accuracy
      EXPECT_EQ(orbits.size(), 31);
      EXPECT_EQ(diff_trans_orbits.size(), 516);
    }
  };

  struct TestOrbits1 : test::TestPrimPeriodicDiffusionTransformationOrbits {
    TestOrbits1(const PrimClex &primclex) :
      test::TestPrimPeriodicDiffusionTransformationOrbits(
        primclex,
        jsonFile(autotools::abs_srcdir() + "/tests/unit/kinetics/ZrO_bspecs_0.json"),
        2, 4) {
      // counts only check for consistency, not accuracy
      EXPECT_EQ(orbits.size(), 74);
      EXPECT_EQ(diff_trans_orbits.size(), 181);
    }
  };

  void
  local_point_clusters_test(const test::TestSupercell &test_supercell,
                            const Kinetics::DiffusionTransformation &phenom) {

    IntegralCluster phenom_sites = phenom.cluster();

    std::vector<PermuteIterator> dt_permute_group =
      phenom.invariant_subgroup(test_supercell.scel);
    SymGroup dt_group = make_sym_group(dt_permute_group, test_supercell.scel.sym_info().supercell_lattice());

    OrbitGenerators<WithinScelIntegralClusterOrbit> generators(
      dt_group, test_supercell.within_scel_sym_compare);

    for(Index l = 0; l < test_supercell.scel.num_sites(); ++l) {
      IntegralCluster el(test_supercell.scel.primclex().prim());
      el.elements().push_back(test_supercell.scel.uccoord(l));
      generators.insert(el);
    }

    std::vector<WithinScelIntegralClusterOrbit> local_point_orbits;
    generators.make_orbits(std::back_inserter(local_point_orbits));

    Index count = 0;
    for(const auto &orbit : local_point_orbits) {
      count += orbit.size();
    }

    EXPECT_EQ(count, test_supercell.scel.num_sites());

    OrbitPrinterOptions opt;
    opt.coord_type = CART;
    FullOrbitPrinter<IntegralCluster> orbit_printer(opt);
    print_clust(local_point_orbits.begin(), local_point_orbits.end(),
                test_supercell.scel.primclex().log(), orbit_printer);
  }

  void
  local_pair_clusters_test(const test::TestSupercell &test_supercell,
                           const Kinetics::DiffusionTransformation &phenom) {

    IntegralCluster phenom_sites = phenom.cluster();

    std::vector<PermuteIterator> dt_permute_group =
      phenom.invariant_subgroup(test_supercell.scel);
    SymGroup dt_group = make_sym_group(dt_permute_group, test_supercell.scel.sym_info().supercell_lattice());

    OrbitGenerators<WithinScelIntegralClusterOrbit> generators(
      dt_group, test_supercell.within_scel_sym_compare);

    for(Index l = 0; l < test_supercell.scel.num_sites(); ++l) {
      for(Index m = 0; m < test_supercell.scel.num_sites(); ++m) {
        if(l < m) {
          IntegralCluster el(test_supercell.scel.primclex().prim());
          el.elements().push_back(test_supercell.scel.uccoord(l));
          el.elements().push_back(test_supercell.scel.uccoord(m));
          generators.insert(el);
        }
      }
    }

    std::vector<WithinScelIntegralClusterOrbit> local_pair_orbits;
    generators.make_orbits(std::back_inserter(local_pair_orbits));

    Index count = 0;
    for(const auto &orbit : local_pair_orbits) {
      count += orbit.size();
    }

    EXPECT_EQ(count, test_supercell.scel.num_sites() *
              (test_supercell.scel.num_sites() - 1) / 2);

    OrbitPrinterOptions opt;
    opt.coord_type = CART;
    FullOrbitPrinter<IntegralCluster> orbit_printer(opt);
    print_clust(local_pair_orbits.begin(), local_pair_orbits.end(),
                test_supercell.scel.primclex().log(), orbit_printer);
  }
}

TEST(ClusterSpecsParserTest, PrimPeriodicClustersByMaxLengthTest) {

  /// Make test project
  test::FCCTernaryProj proj;
  proj.check_init();
  PrimClex primclex(proj.dir, null_log());
  EXPECT_EQ(true, true);

  auto &log = null_log();
  {
    auto parser = PrimParser(primclex, R"({
      "orbit_branch_specs": {
        "2": {"max_length": 3.}
      }
    })");

    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->orbit_branch_specs().warning.size(), 0);

    EXPECT_EQ(parser->max_branch(), 2);
    EXPECT_EQ(almost_equal(parser->max_length(2), 3.0), true);
    EXPECT_THROW(parser->max_length(1), std::invalid_argument);
    EXPECT_THROW(parser->max_length(3), std::invalid_argument);
    EXPECT_EQ(parser->custom_generators().elements.size(), 0);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  {
    auto parser = PrimParser(primclex, R"({})");

    EXPECT_EQ(parser->valid(), false);
    EXPECT_EQ(parser->error.size(), 1);
    EXPECT_EQ(parser->warning.size(), 0);
    EXPECT_EQ(parser->orbit_branch_specs().warning.size(), 0);

    EXPECT_EQ(parser->max_branch(), 1);
    EXPECT_EQ(parser->custom_generators().elements.size(), 0);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  {
    auto parser = PrimParser(primclex, R"({"orbit_branch_specs": []})");

    EXPECT_EQ(parser->valid(), false);
    EXPECT_EQ(parser->orbit_branch_specs().error.size(), 1);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  {
    auto parser = PrimParser(primclex, R"({
      "orbit_branch_specs": {
        "1": {"max_length": 3.},
        "2": {"max_length": 3.}
      }
    })");

    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->orbit_branch_specs().warning.size(), 1);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  {
    auto parser = PrimParser(primclex, R"({
      "orbit_branch_specs": {
        "3": {"max_length": 3.}
      }
    })");

    EXPECT_EQ(parser->valid(), false);
    EXPECT_EQ(parser->orbit_branch_specs().warning.size(), 0);
    EXPECT_EQ(parser->orbit_branch_specs().error.size(), 1);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  {
    auto parser = PrimParser(primclex, R"({
      "orbit_branch_specs": {
        "2": {"max_length": 3.},
        "3": {"max_length": 4.}
      }
    })");

    EXPECT_EQ(parser->valid(), false);
    EXPECT_EQ(parser->orbit_branch_specs().warning.size(), 0);
    EXPECT_EQ(parser->orbit_branch_specs().error.size(), 1);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  {
    auto parser = PrimParser(primclex, R"({
      "orbit_branch_specs": {
        "2": {"max_length": 4.},
        "3": {"max_length": 4.}
      }
    })");

    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->orbit_branch_specs().warning.size(), 0);
    EXPECT_EQ(parser->orbit_branch_specs().error.size(), 0);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  {
    auto parser = PrimParser(primclex, R"({
      "orbit_branch_specs": {
        "2": {"max_length": 4., "something_else": 1},
        "3": {"max_length": 4.},
        "_3": {"max_length": 5.}
      }
    })");

    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->orbit_branch_specs().warning.size(), 2);
    EXPECT_EQ(parser->orbit_branch_specs().error.size(), 0);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  // --- PrimPeriodicOrbitSpecsParser ----

  // --- valid, w/ subclusters ---
  {
    // includes null cluster, point cluster, pair cluster == 3 custom generators
    auto parser = PrimParser(primclex, R"({
      "orbit_specs": [
        {
          "coordinate_mode": "Integral",
          "prototype": [
            [0, 0, 0, 0],
            [0, 1, 0, 0]
          ]
        }
      ]
    })");

    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->orbit_specs().warning.size(), 0);
    EXPECT_EQ(parser->orbit_specs().error.size(), 0);
    EXPECT_EQ(parser->custom_generators().elements.size(), 3);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  // --- valid, no subclusters ---
  {
    auto parser = PrimParser(primclex, R"({
      "orbit_specs": [
        {
          "coordinate_mode": "Integral",
          "prototype": [
            [0, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0]
          ],
          "include_subclusters": false
        }
      ]
    })");

    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->orbit_specs().warning.size(), 0);
    EXPECT_EQ(parser->orbit_specs().error.size(), 0);
    EXPECT_EQ(parser->custom_generators().elements.size(), 1);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  // --- valid, include_subclusters ---
  {
    auto parser = PrimParser(primclex, R"({
      "orbit_specs": [
        {
          "coordinate_mode": "Integral",
          "prototype": [
            [0, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0]
          ],
          "include_subclusters": true
        }
      ]
    })");

    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->orbit_specs().warning.size(), 0);
    EXPECT_EQ(parser->orbit_specs().error.size(), 0);
    EXPECT_EQ(parser->custom_generators().elements.size(), 4);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  // --- valid w/warning ---
  {
    auto parser = PrimParser(primclex, R"({
      "orbit_specs": [
        {
          "coordinate_mode": "Integral",
          "prototype": [
            [0, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0]
          ],
          "something_else": 1
        }
      ]
    })");

    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->orbit_specs().warning.size(), 1);
    EXPECT_EQ(parser->orbit_specs().error.size(), 0);
    EXPECT_EQ(parser->custom_generators().elements.size(), 4);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  // --- error, not array ---
  {
    auto parser = PrimParser(primclex, R"({
      "orbit_specs": {
        "coordinate_mode": "Integral",
        "prototype": [
          [0, 0, 0, 0],
          [0, 1, 0, 0]
        ]
      }
    })");

    EXPECT_EQ(parser->valid(), false);
    EXPECT_EQ(parser->orbit_specs().error.size(), 1);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }
}

TEST(ClusterSpecsParserTest, LocalClustersByMaxLengthTest) {
  /// Make test project
  test::FCCTernaryProj proj;
  proj.check_init();
  Log &log = null_log();
  PrimClex primclex(proj.dir, null_log());
  EXPECT_EQ(true, true);

  TestOrbits0 to(primclex);

  jsonFile diff_trans_json {autotools::abs_srcdir() + "/tests/unit/kinetics/FCCTernary_diff_trans_0.json"};
  Completer::EnumOption enum_opt;
  enum_opt.desc();
  int success = Kinetics::DiffusionTransformationEnum::run(primclex, diff_trans_json, enum_opt, nullptr);
  const auto &diff_trans_db = primclex.generic_db<Kinetics::PrimPeriodicDiffTransOrbit>();
  EXPECT_EQ(diff_trans_db.size(), 28);
  EXPECT_EQ(success, 0);

  OrbitPrinterOptions opt;
  opt.coord_type = CART;
  PrototypePrinter<Kinetics::DiffusionTransformation> dt_orbit_printer(opt);

  //  {
  //    Index i = 0;
  //    for(const auto &orbit : diff_trans_db) {
  //      log << "\n --- " << orbit.name() << " ---" << std::endl;
  //      dt_orbit_printer(orbit, log, i++, diff_trans_db.size());
  //    }
  //  }

  // generate local orbits (w/ scel symmetry) using a 3x3x3 standard FCC supercell (so 3x3x3x4 sites)
  test::TestStandardFCCSupercell fcc_supercell(primclex, 3, 3, 3);

  //print prototypes
  //  {
  //    PrototypePrinter<IntegralCluster> printer;
  //    print_clust(to.orbits.begin(), to.orbits.end(), log, printer);
  //  }

  //print DiffTrans prototypes
  //  {
  //    PrototypePrinter<Kinetics::DiffusionTransformation> printer;
  //    print_clust(to.diff_trans_orbits.begin(), to.diff_trans_orbits.end(), log, printer);
  //    jsonParser json;
  //    write_clust(to.diff_trans_orbits.begin(), to.diff_trans_orbits.end(), json, printer);
  //    log << json << std::endl;
  //  }

  // --- minimal ---
  {
    auto parser = LocalParser(fcc_supercell.scel, R"({
      "standard": {
        "orbit_branch_specs": {
          "1": {"cutoff_radius": 10.}
        }
      }
    })");

    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->standard->max_branch, 1);
    EXPECT_EQ(parser->standard->max_length_including_phenomenal(), false);
    EXPECT_EQ(almost_equal(parser->standard->cutoff_radius(1), 10.), true);
    EXPECT_THROW(parser->standard->max_length(1), std::runtime_error);

    EXPECT_EQ(parser->all_warnings().size(), 0);
    EXPECT_EQ(parser->all_errors().size(), 0);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  // --- typical ---
  {
    auto parser = LocalParser(fcc_supercell.scel, R"({
      "standard": {
        "orbit_branch_specs": {
          "1": {"cutoff_radius": 3.},
          "2": {"max_length": 12.0, "cutoff_radius": 3.},
          "3": {"max_length": 10.0, "cutoff_radius": 3.}
        }
      }
    })");

    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->standard->max_branch, 3);
    EXPECT_EQ(parser->standard->max_length_including_phenomenal(), false);
    EXPECT_EQ(almost_equal(parser->standard->cutoff_radius(1), 3.), true);
    EXPECT_EQ(almost_equal(parser->standard->cutoff_radius(2), 3.), true);
    EXPECT_EQ(almost_equal(parser->standard->cutoff_radius(3), 3.), true);
    EXPECT_THROW(parser->standard->max_length(1), std::runtime_error);
    EXPECT_EQ(almost_equal(parser->standard->max_length(2), 12.), true);
    EXPECT_EQ(almost_equal(parser->standard->max_length(3), 10.), true);

    EXPECT_EQ(parser->all_warnings().size(), 0);
    EXPECT_EQ(parser->all_errors().size(), 0);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  // --- typical, warning for max_length on brach 1 w/out max_length_including_phenomenal option ---
  {
    auto parser = LocalParser(fcc_supercell.scel, R"({
      "standard": {
        "orbit_branch_specs": {
          "1": {"max_length": 12.0, "cutoff_radius": 3.},
          "2": {"max_length": 12.0, "cutoff_radius": 3.},
          "3": {"max_length": 10.0, "cutoff_radius": 3.}
        }
      }
    })");

    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->standard->warning.size(), 1);
    EXPECT_EQ(parser->standard->error.size(), 0);
    EXPECT_EQ(parser->all_warnings().size(), 1);
    EXPECT_EQ(parser->all_errors().size(), 0);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  // --- typical, error for nonincreasing ---
  {
    auto parser = LocalParser(fcc_supercell.scel, R"({
      "standard": {
        "orbit_branch_specs": {
          "1": {"cutoff_radius": 3.},
          "2": {"max_length": 12.0, "cutoff_radius": 4.},
          "3": {"max_length": 13.0, "cutoff_radius": 3.}
        }
      }
    })");

    EXPECT_EQ(parser->valid(), false);
    EXPECT_EQ(parser->standard->warning.size(), 0);
    EXPECT_EQ(parser->standard->error.size(), 2);
    EXPECT_EQ(parser->all_warnings().size(), 0);
    EXPECT_EQ(parser->all_errors().size(), 2);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  // --- typical, w/ max_length_including_phenomenal ---
  {
    auto parser = LocalParser(fcc_supercell.scel, R"({
      "standard": {
        "orbit_branch_specs": {
          "max_length_including_phenomenal": true,
          "1": {"max_length": 12.0, "cutoff_radius": 3.},
          "2": {"max_length": 12.0, "cutoff_radius": 3.},
          "3": {"max_length": 10.0, "cutoff_radius": 3.}
        }
      }
    })");

    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->standard->max_branch, 3);
    EXPECT_EQ(parser->standard->max_length_including_phenomenal(), true);
    EXPECT_EQ(almost_equal(parser->standard->cutoff_radius(1), 3.), true);
    EXPECT_EQ(almost_equal(parser->standard->cutoff_radius(2), 3.), true);
    EXPECT_EQ(almost_equal(parser->standard->cutoff_radius(3), 3.), true);
    EXPECT_EQ(almost_equal(parser->standard->max_length(1), 12.), true);
    EXPECT_EQ(almost_equal(parser->standard->max_length(2), 12.), true);
    EXPECT_EQ(almost_equal(parser->standard->max_length(3), 10.), true);

    EXPECT_EQ(parser->all_warnings().size(), 0);
    EXPECT_EQ(parser->all_errors().size(), 0);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  // --- typical, w/ max_length_including_phenomenal && max_length="inf" (all clusters of branch in neighborhood) ---
  {
    auto parser = LocalParser(fcc_supercell.scel, R"({
      "standard": {
        "orbit_branch_specs": {
          "max_length_including_phenomenal": true,
          "1": {"max_length": "inf", "cutoff_radius": 3.},
          "2": {"max_length": "inf", "cutoff_radius": 3.},
          "3": {"max_length": "inf", "cutoff_radius": 3.}
        }
      }
    })");

    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->standard->max_branch, 3);
    EXPECT_EQ(parser->standard->max_length_including_phenomenal(), true);
    EXPECT_EQ(almost_equal(parser->standard->cutoff_radius(1), 3.), true);
    EXPECT_EQ(almost_equal(parser->standard->cutoff_radius(2), 3.), true);
    EXPECT_EQ(almost_equal(parser->standard->cutoff_radius(3), 3.), true);
    EXPECT_EQ(parser->standard->max_length(1), std::numeric_limits<double>::infinity());
    EXPECT_EQ((3.0 < parser->standard->max_length(1)), true);
    EXPECT_EQ(parser->standard->max_length(2), std::numeric_limits<double>::infinity());
    EXPECT_EQ(parser->standard->max_length(3), std::numeric_limits<double>::infinity());

    EXPECT_EQ(parser->all_warnings().size(), 0);
    EXPECT_EQ(parser->all_errors().size(), 0);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  // --- w/ max_length_including_phenomenal, error for nonincreasing ---
  {
    auto parser = LocalParser(fcc_supercell.scel, R"({
      "standard": {
        "orbit_branch_specs": {
          "max_length_including_phenomenal": true,
          "1": {"max_length": 13.0, "cutoff_radius": 3.},
          "2": {"max_length": 12.0, "cutoff_radius": 4.},
          "3": {"max_length": 13.0, "cutoff_radius": 3.}
        }
      }
    })");

    EXPECT_EQ(parser->valid(), false);
    EXPECT_EQ(parser->standard->warning.size(), 0);
    EXPECT_EQ(parser->standard->error.size(), 2);
    EXPECT_EQ(parser->all_warnings().size(), 0);
    EXPECT_EQ(parser->all_errors().size(), 2);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  // --- custom orbit_branch_specs, equivalence_type=prim (default)  ---
  {
    auto parser = LocalParser(fcc_supercell.scel, R"({
      "standard": {
        "orbit_branch_specs": {
          "max_length_including_phenomenal": true,
          "1": {"max_length": 4., "cutoff_radius": 3.},
          "2": {"max_length": 4., "cutoff_radius": 3.}
        }
      },
      "custom": [
        {
          "phenomenal": {
            "occ_transform" : [
              {"from_value" : 2, "to_value" : 2, "uccoord" : [ 0, 0, 0, 0 ]},
              {"from_value" : 2, "to_value" : 2, "uccoord" : [ 0, 0, 0, 1 ]}
            ],
            "species_traj" : [
              {
                "from" : {"occ" : 2, "pos" : 0, "uccoord" : [ 0, 0, 0, 0 ]},
                "to" : {"occ" : 2, "pos" : 0, "uccoord" : [ 0, 0, 0, 1 ]}
              },
              {
                "from" : {"occ" : 2, "pos" : 0, "uccoord" : [ 0, 0, 0, 1 ]},
                "to" : {"occ" : 2, "pos" : 0, "uccoord" : [ 0, 0, 0, 0 ]}
              }
            ]
          },
          "orbit_branch_specs": {
            "max_length_including_phenomenal": true,
            "1": {"max_length": "inf", "cutoff_radius": 6.},
            "2": {"max_length": "inf", "cutoff_radius": 6.},
            "3": {"max_length": "inf", "cutoff_radius": 6.}
          }
        }
      ]
    })");

    EXPECT_EQ(parser->valid(), true);

    EXPECT_EQ(parser->standard->max_branch, 2);
    EXPECT_EQ(parser->standard->max_length_including_phenomenal(), true);
    EXPECT_EQ(almost_equal(parser->standard->cutoff_radius(1), 3.), true);
    EXPECT_EQ(almost_equal(parser->standard->cutoff_radius(2), 3.), true);
    EXPECT_EQ(almost_equal(parser->standard->max_length(1), 4.), true);
    EXPECT_EQ(almost_equal(parser->standard->max_length(2), 4.), true);

    EXPECT_EQ(parser->custom->data.size(), 1);
    EXPECT_EQ(parser->custom->data[0].orbit_branch_specs->max_branch, 3);
    EXPECT_EQ(parser->custom->data[0].orbit_branch_specs->max_length_including_phenomenal(), true);
    EXPECT_EQ(almost_equal(parser->custom->data[0].orbit_branch_specs->cutoff_radius(1), 6.), true);
    EXPECT_EQ(almost_equal(parser->custom->data[0].orbit_branch_specs->cutoff_radius(2), 6.), true);
    EXPECT_EQ(almost_equal(parser->custom->data[0].orbit_branch_specs->cutoff_radius(3), 6.), true);
    EXPECT_EQ(parser->custom->data[0].orbit_branch_specs->max_length(1), std::numeric_limits<double>::infinity());
    EXPECT_EQ(parser->custom->data[0].orbit_branch_specs->max_length(2), std::numeric_limits<double>::infinity());
    EXPECT_EQ(parser->custom->data[0].orbit_branch_specs->max_length(3), std::numeric_limits<double>::infinity());

    // checks comparing input phenomenal cluster to test cluster
    // the input custom phenomenal cluster is to.diff_trans_orbits[0].prototype(), so
    // only clusters in orbit to.diff_trans_orbits[0] should be found
    if(parser->valid()) {
      Index linear_orbit_index = 0;
      for(const auto &orbit : to.diff_trans_orbits) {
        for(const auto &equiv : orbit) {
          auto find_res = parser->find(equiv);
          auto it = find_res.first;
          auto op = find_res.second;
          EXPECT_EQ((it != parser->custom->data.end()), (linear_orbit_index == 0));

          // check op*phenom == equiv
          if(it != parser->custom->data.end()) {
            auto &f = *it->phenom->prim_sym_compare;
            EXPECT_EQ(f.equal(f.prepare(copy_apply(op, *it->phenom->phenom)), f.prepare(equiv)), true);
          }

          int expected_max_branch = (linear_orbit_index == 0 ? 3 : 2);
          EXPECT_EQ(parser->max_branch(it), expected_max_branch);
          EXPECT_EQ(parser->max_length_including_phenomenal(it), true);

          for(int branch = 1; branch <= parser->max_branch(it); ++branch) {
            double expected_cutoff_radius = (linear_orbit_index == 0 ? 6. : 3.);
            EXPECT_EQ(parser->cutoff_radius(it, branch), expected_cutoff_radius);

            double expected_max_length = (linear_orbit_index == 0 ? std::numeric_limits<double>::infinity() : 4.);
            EXPECT_EQ(parser->max_length(it, branch), expected_max_length);
          }
        }
        linear_orbit_index++;
      }
    }

    EXPECT_EQ(parser->all_warnings().size(), 0);
    EXPECT_EQ(parser->all_errors().size(), 0);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

  // --- custom orbit_branch_specs, equivalence_type=prim (default), specified by name  ---
  {
    //    auto parser = LocalParser(ts.scel, R"({
    //      "standard": {
    //        "orbit_branch_specs": {
    //          "max_length_including_phenomenal": true,
    //          "1": {"max_length": 4., "cutoff_radius": 3.},
    //          "2": {"max_length": 4., "cutoff_radius": 3.}
    //        }
    //      },
    //      "custom": [
    //        {
    //          "phenomenal": "<name>",
    //          "orbit_branch_specs": {
    //            "max_length_including_phenomenal": true,
    //            "1": {"max_length": "inf", "cutoff_radius": 6.},
    //            "2": {"max_length": "inf", "cutoff_radius": 6.},
    //            "3": {"max_length": "inf", "cutoff_radius": 6.}
    //          }
    //        }
    //      ]
    //    })");
    //
    //    EXPECT_EQ(parser->valid(), true);
    //    EXPECT_EQ(parser->standard->max_branch, 2);
    //    EXPECT_EQ(parser->standard->max_length_including_phenomenal(), true);
    //    EXPECT_EQ(almost_equal(parser->standard->cutoff_radius(1), 3.), true);
    //    EXPECT_EQ(almost_equal(parser->standard->cutoff_radius(2), 3.), true);
    //    EXPECT_EQ(almost_equal(parser->standard->max_length(1), 4.), true);
    //    EXPECT_EQ(almost_equal(parser->standard->max_length(2), 4.), true);
    //
    //    EXPECT_EQ(parser->custom->data.size(), 1);
    //    EXPECT_EQ(parser->custom->data[0].orbit_branch_specs->max_branch, 3);
    //    EXPECT_EQ(parser->custom->data[0].orbit_branch_specs->max_length_including_phenomenal(), true);
    //    EXPECT_EQ(almost_equal(parser->custom->data[0].orbit_branch_specs->cutoff_radius(1), 6.), true);
    //    EXPECT_EQ(almost_equal(parser->custom->data[0].orbit_branch_specs->cutoff_radius(2), 6.), true);
    //    EXPECT_EQ(almost_equal(parser->custom->data[0].orbit_branch_specs->cutoff_radius(3), 6.), true);
    //    EXPECT_EQ(parser->custom->data[0].orbit_branch_specs->max_length(1), std::numeric_limits<double>::infinity());
    //    EXPECT_EQ(parser->custom->data[0].orbit_branch_specs->max_length(2), std::numeric_limits<double>::infinity());
    //    EXPECT_EQ(parser->custom->data[0].orbit_branch_specs->max_length(3), std::numeric_limits<double>::infinity());
    //
    //    // checks comparing input phenomenal cluster to test cluster
    //    // the input custom phenomenal cluster is to.diff_trans_orbits[0].prototype(), so
    //    // only clusters in orbit to.diff_trans_orbits[0] should be found
    //    {
    //      Index linear_orbit_index = 0;
    //      for(const auto &orbit : to.diff_trans_orbits) {
    //        for(const auto &equiv : orbit) {
    //          auto find_res = parser->find(equiv);
    //          auto it = find_res.first;
    //          auto op = find_res.second;
    //          EXPECT_EQ((it != parser->custom->data.end()), (linear_orbit_index == 0));
    //
    //          // check op*phenom == equiv
    //          if(it != parser->custom->data.end()) {
    //            auto &f = *it->phenom->prim_sym_compare;
    //            EXPECT_EQ(f.equal(f.prepare(copy_apply(op, *it->phenom->phenom)), f.prepare(equiv)), true);
    //          }
    //
    //          int expected_max_branch = (linear_orbit_index == 0 ? 3 : 2);
    //          EXPECT_EQ(parser->max_branch(it), expected_max_branch);
    //          EXPECT_EQ(parser->max_length_including_phenomenal(it), true);
    //
    //          for(int branch = 1; branch <= parser->max_branch(it); ++branch) {
    //            double expected_cutoff_radius = (linear_orbit_index == 0 ? 6. : 3.);
    //            EXPECT_EQ(parser->cutoff_radius(it, branch), expected_cutoff_radius);
    //
    //            double expected_max_length = (linear_orbit_index == 0 ? std::numeric_limits<double>::infinity() : 4.);
    //            EXPECT_EQ(parser->max_length(it, branch), expected_max_length);
    //          }
    //        }
    //        linear_orbit_index++;
    //      }
    //    }
    //
    //    EXPECT_EQ(parser->all_warning().size(), 0);
    //    EXPECT_EQ(parser->all_error().size(), 0);
    //
    //    log << parser->report() << std::endl;
  }

  // --- custom orbit_specs, equivalence_type=prim (default)  ---
  {
    auto parser = LocalParser(fcc_supercell.scel, R"({
      "standard": {
        "orbit_branch_specs": {
          "max_length_including_phenomenal": true,
          "1": {"max_length": 4., "cutoff_radius": 3.},
          "2": {"max_length": 4., "cutoff_radius": 3.}
        }
      },
      "custom": [
        {
          "phenomenal": {
            "occ_transform" : [
              {"from_value" : 2, "to_value" : 2, "uccoord" : [ 0, 0, 0, 0 ]},
              {"from_value" : 2, "to_value" : 2, "uccoord" : [ 0, 0, 0, 1 ]}
            ],
            "species_traj" : [
              {
                "from" : {"occ" : 2, "pos" : 0, "uccoord" : [ 0, 0, 0, 0 ]},
                "to" : {"occ" : 2, "pos" : 0, "uccoord" : [ 0, 0, 0, 1 ]}
              },
              {
                "from" : {"occ" : 2, "pos" : 0, "uccoord" : [ 0, 0, 0, 1 ]},
                "to" : {"occ" : 2, "pos" : 0, "uccoord" : [ 0, 0, 0, 0 ]}
              }
            ]
          },
          "orbit_specs": [
            {
              "coordinate_mode": "Integral",
              "prototype": [
                [0, 1, 0, 0]
              ]
            }
          ]
        }
      ]
    })");

    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->standard->max_branch, 2);
    EXPECT_EQ(parser->standard->max_length_including_phenomenal(), true);
    EXPECT_EQ(almost_equal(parser->standard->cutoff_radius(1), 3.), true);
    EXPECT_EQ(almost_equal(parser->standard->cutoff_radius(2), 3.), true);
    EXPECT_EQ(almost_equal(parser->standard->max_length(1), 4.), true);
    EXPECT_EQ(almost_equal(parser->standard->max_length(2), 4.), true);

    EXPECT_EQ(parser->custom->data.size(), 1);
    EXPECT_EQ(parser->custom->data[0].orbit_specs->prototypes.size(), 1);
    EXPECT_EQ(parser->custom->data[0].orbit_specs->prototypes[0].include_subclusters, 1);
    EXPECT_EQ(parser->custom->data[0].orbit_specs->prototypes[0].cluster.size(), 1);

    ScelPeriodicDiffTransSymCompare dt_scel_sym_compare(
      fcc_supercell.scel.primclex().shared_prim(),
      fcc_supercell.scel.transf_mat(),
      fcc_supercell.scel.crystallography_tol());
    // checks comparing input phenomenal cluster to test cluster
    // the input custom phenomenal cluster is to.diff_trans_orbits[0].prototype(), so
    // only clusters in orbit to.diff_trans_orbits[0] should be found
    if(parser->valid()) {
      Index linear_orbit_index = 0;
      for(const auto &orbit : to.diff_trans_orbits) {
        Index equiv_index = 0;
        for(const auto &equiv : orbit) {

          auto find_res = parser->find(equiv);
          auto it = find_res.first;
          auto op = find_res.second;

          // only linear_orbit_index == 0 should match custom orbit specs
          EXPECT_EQ((it != parser->custom->data.end()), (linear_orbit_index == 0));

          // generate local orbit
          if(it != parser->custom->data.end()) {

            // check op*phenom == equiv
            auto &f = *it->phenom->prim_sym_compare;
            EXPECT_EQ(f.equal(f.prepare(copy_apply(op, *it->phenom->phenom)), f.prepare(equiv)), true);

            std::vector<PermuteIterator> dt_permute_group =
              equiv.invariant_subgroup(fcc_supercell.scel);
            SymGroup dt_group = make_sym_group(dt_permute_group, fcc_supercell.scel.sym_info().supercell_lattice());
            OrbitGenerators<ScelPeriodicIntegralClusterOrbit> test_gen(
              dt_group, fcc_supercell.scel_sym_compare);
            parser->insert_custom_generators(find_res, test_gen, sym::CopyApplyWithPrim_f(primclex.shared_prim()));

            if(linear_orbit_index == 0) {

              EXPECT_EQ(dt_group.size(), 8);

              // 2 elements (because default==include_subclusters -> null cluster)
              EXPECT_EQ(test_gen.elements.size(), 2);

              auto gen_it = test_gen.elements.begin();
              {
                // orbit of null cluster -> 1 equiv
                ScelPeriodicIntegralClusterOrbit orbit(
                  *gen_it++, dt_group, fcc_supercell.scel_sym_compare);
                EXPECT_EQ(orbit.size(), 1);
              }

              {
                // should generate 4 equivalents (there are 4NN of both sites in 1NN pair cluster)
                ScelPeriodicIntegralClusterOrbit orbit(
                  *gen_it++, dt_group, fcc_supercell.scel_sym_compare);
                EXPECT_EQ(orbit.size(), 4);
              }
            }
            else {
              // pass
            }
          }
          equiv_index++;
        }
        linear_orbit_index++;
      }
    }

    EXPECT_EQ(parser->all_warnings().size(), 0);
    EXPECT_EQ(parser->all_errors().size(), 0);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }
}

TEST(ClusterSpecsParserTest, LocalClustersByMaxLengthTest_Tet) {
  /// Make test project
  test::FCCTernaryProj proj;
  proj.check_init();
  Log &log = null_log();
  PrimClex primclex(proj.dir, null_log());
  EXPECT_EQ(true, true);

  TestOrbits0 to(primclex);


  // generate local orbits (w/ scel symmetry) using a 4x3x3 standard FCC supercell (so 4x3x3x4 sites)
  test::TestStandardFCCSupercell fcc_supercell(primclex, 4, 3, 3);

  // --- custom orbit_specs, equivalence_type=scel, non-cubic supercell  ---
  {
    auto parser = LocalParser(fcc_supercell.scel, R"({
      "standard": {
        "orbit_branch_specs": {
          "max_length_including_phenomenal": true,
          "1": {"max_length": 4., "cutoff_radius": 3.},
          "2": {"max_length": 4., "cutoff_radius": 3.}
        }
      },
      "custom": [
        {
          "phenomenal": {
            "equivalence_type": "scel",
            "occ_transform" : [
              {"from_value" : 2, "to_value" : 2, "uccoord" : [ 0, 0, 0, 0 ]},
              {"from_value" : 2, "to_value" : 2, "uccoord" : [ 0, 0, 0, 1 ]}
            ],
            "species_traj" : [
              {
                "from" : {"occ" : 2, "pos" : 0, "uccoord" : [ 0, 0, 0, 0 ]},
                "to" : {"occ" : 2, "pos" : 0, "uccoord" : [ 0, 0, 0, 1 ]}
              },
              {
                "from" : {"occ" : 2, "pos" : 0, "uccoord" : [ 0, 0, 0, 1 ]},
                "to" : {"occ" : 2, "pos" : 0, "uccoord" : [ 0, 0, 0, 0 ]}
              }
            ]
          },
          "orbit_specs": [
            {
              "coordinate_mode": "Integral",
              "prototype": [
                [0, 1, 0, 0]
              ]
            }
          ]
        },
        {
          "phenomenal": {
            "equivalence_type": "scel",
            "occ_transform" : [
              {"from_value" : 2, "to_value" : 2, "uccoord" : [ 0, 0, 0, 0 ]},
              {"from_value" : 2, "to_value" : 2, "uccoord" : [ 0, 1, 0, 0 ]}
            ],
            "species_traj" : [
              {
                "from" : {"occ" : 2, "pos" : 0, "uccoord" : [ 0, 0, 0, 0 ]},
                "to" : {"occ" : 2, "pos" : 0, "uccoord" : [ 0, 1, 0, 0 ]}
              },
              {
                "from" : {"occ" : 2, "pos" : 0, "uccoord" : [ 0, 1, 0, 0 ]},
                "to" : {"occ" : 2, "pos" : 0, "uccoord" : [ 0, 0, 0, 0 ]}
              }
            ]
          },
          "orbit_specs": [
            {
              "coordinate_mode": "Integral",
              "prototype": [
                [0, 0, 0, 1]
              ]
            }
          ]
        }
      ]
    })");

    EXPECT_EQ(parser->valid(), true);
    EXPECT_EQ(parser->standard->max_branch, 2);
    EXPECT_EQ(parser->standard->max_length_including_phenomenal(), true);
    EXPECT_EQ(almost_equal(parser->standard->cutoff_radius(1), 3.), true);
    EXPECT_EQ(almost_equal(parser->standard->cutoff_radius(2), 3.), true);
    EXPECT_EQ(almost_equal(parser->standard->max_length(1), 4.), true);
    EXPECT_EQ(almost_equal(parser->standard->max_length(2), 4.), true);

    EXPECT_EQ(parser->custom->data.size(), 2);
    EXPECT_EQ(parser->custom->data[0].orbit_specs->prototypes.size(), 1);
    EXPECT_EQ(parser->custom->data[0].orbit_specs->prototypes[0].include_subclusters, 1);
    EXPECT_EQ(parser->custom->data[0].orbit_specs->prototypes[0].cluster.size(), 1);
    EXPECT_EQ(parser->custom->data[1].orbit_specs->prototypes.size(), 1);
    EXPECT_EQ(parser->custom->data[1].orbit_specs->prototypes[0].include_subclusters, 1);
    EXPECT_EQ(parser->custom->data[1].orbit_specs->prototypes[0].cluster.size(), 1);

    ScelPeriodicDiffTransSymCompare diff_trans_scel_sym_compare(
      fcc_supercell.scel.primclex().shared_prim(),
      fcc_supercell.scel.transf_mat(),
      fcc_supercell.scel.crystallography_tol());
    // checks comparing input phenomenal cluster to test cluster
    // the input custom phenomenal cluster is to.diff_trans_orbits[0].prototype(), so
    // 4/6 clusters in orbit to.diff_trans_orbits[0] should match the first custom phenom
    // 2/6 should match the other custom phenom

    OrbitPrinterOptions opt;
    opt.coord_type = CART;
    Printer<IntegralCluster> cluster_printer(opt);

    //    log << "check # matches" << std::endl;
    if(parser->valid()) {
      std::vector<Index> count = {0, 0};
      EXPECT_EQ(to.diff_trans_orbits[0].size(), 6);
      for(const auto &equiv : to.diff_trans_orbits[0]) {

        //        log << "\n##################" << std::endl;
        //        log << "equiv: \n" << equiv << std::endl;
        //        log << "equiv.cluster(): \n";
        cluster_printer.print(equiv.cluster(), log);
        log << std::endl;

        auto find_res = parser->find(equiv);
        auto it = find_res.first;
        auto op = find_res.second;
        if(it != parser->custom->data.end()) {
          auto &f = *it->phenom->scel_sym_compare;
          EXPECT_EQ(f.equal(f.prepare(copy_apply(op, *it->phenom->phenom)), f.prepare(equiv)), true);

          std::vector<PermuteIterator> dt_permute_group =
            equiv.invariant_subgroup(fcc_supercell.scel);
          SymGroup dt_group = make_sym_group(dt_permute_group, fcc_supercell.scel.sym_info().supercell_lattice());

          if(&*it == &parser->custom->data[0]) {
            count[0]++;
            EXPECT_EQ(dt_group.size(), 4);
          }
          else if(&*it == &parser->custom->data[1]) {
            count[1]++;
            EXPECT_EQ(dt_group.size(), 8);
          }

          OrbitGenerators<ScelPeriodicIntegralClusterOrbit> test_gen(
            dt_group, fcc_supercell.scel_sym_compare);
          parser->insert_custom_generators(find_res, test_gen, sym::CopyApplyWithPrim_f(primclex.shared_prim()));



          // 2 elements (because default==include_subclusters -> null cluster)
          EXPECT_EQ(test_gen.elements.size(), 2);

          auto gen_it = test_gen.elements.begin();
          {
            // orbit of null cluster -> 1 equiv
            ScelPeriodicIntegralClusterOrbit orbit(
              *gen_it++, dt_group, fcc_supercell.scel_sym_compare);
            EXPECT_EQ(orbit.size(), 1);
          }

          {
            // should generate 4 equivalents (the 4 NN of both sites in 1NN pair cluster)
            ScelPeriodicIntegralClusterOrbit orbit(
              *gen_it++, dt_group, fcc_supercell.scel_sym_compare);
            EXPECT_EQ(orbit.size(), 4);
          }
        }
      }
      EXPECT_EQ(count[0], 4); // 4/6 equiv match first phenom
      EXPECT_EQ(count[1], 2); // 2/6 equiv match second phenom
    }

    // other orbits should not match (this is slow)
    log << "check # not matches" << std::endl;
    {
      for(Index i = 1; i < to.diff_trans_orbits.size(); ++i) {
        for(const auto &equiv : to.diff_trans_orbits[i]) {
          auto find_res = parser->find(equiv);
          auto it = find_res.first;
          auto op = find_res.second;
          EXPECT_EQ((it == parser->custom->data.end()), true);
        }
      }
    }

    EXPECT_EQ(parser->all_warnings().size(), 0);
    EXPECT_EQ(parser->all_errors().size(), 0);

    log.begin("Report");
    log << parser->report() << std::endl;
    log << std::endl;
  }

}

TEST(ClusterSpecsParserTest, LocalPointClustersByMaxLengthTest_FCCTernary) {
  /// Make test project
  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  EXPECT_EQ(true, true);

  TestOrbits0 to(primclex);
  const auto &phenom = to.diff_trans_orbits[0].prototype();

  {
    test::TestStandardFCCSupercell ts(primclex, 3, 3, 3);
    local_point_clusters_test(ts, phenom);
  }

  {
    test::TestStandardFCCSupercell ts(primclex, 4, 3, 3);
    local_point_clusters_test(ts, phenom);
  }

  {
    test::TestStandardFCCSupercell ts(primclex, 5, 4, 3);
    local_point_clusters_test(ts, phenom);
  }


  {
    test::TestStandardFCCSupercell ts(primclex, 1, 1, 1);
    local_pair_clusters_test(ts, phenom);
  }

  {
    test::TestStandardFCCSupercell ts(primclex, 2, 2, 2);
    local_pair_clusters_test(ts, phenom);
  }

  {
    test::TestStandardFCCSupercell ts(primclex, 3, 2, 2);
    local_pair_clusters_test(ts, phenom);
  }

  {
    test::TestStandardFCCSupercell ts(primclex, 4, 3, 2);
    local_pair_clusters_test(ts, phenom);
  }
}

TEST(ClusterSpecsParserTest, LocalPointClustersByMaxLengthTest_ZrO) {
  /// Make test project
  test::ZrOProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  EXPECT_EQ(true, true);

  TestOrbits1 to(primclex);
  const auto &phenom = to.diff_trans_orbits[0].prototype();

  {
    test::TestSupercell ts(primclex, Eigen::Vector3i(3, 3, 3).asDiagonal());
    local_point_clusters_test(ts, phenom);
  }

  {
    test::TestSupercell ts(primclex, Eigen::Vector3i(4, 3, 3).asDiagonal());
    local_point_clusters_test(ts, phenom);
  }

  {
    test::TestSupercell ts(primclex, Eigen::Vector3i(5, 4, 3).asDiagonal());
    local_point_clusters_test(ts, phenom);
  }


  {
    test::TestSupercell ts(primclex, Eigen::Vector3i(1, 1, 1).asDiagonal());
    local_pair_clusters_test(ts, phenom);
  }

  {
    test::TestSupercell ts(primclex, Eigen::Vector3i(2, 2, 2).asDiagonal());
    local_pair_clusters_test(ts, phenom);
  }

  {
    test::TestSupercell ts(primclex, Eigen::Vector3i(3, 2, 2).asDiagonal());
    local_pair_clusters_test(ts, phenom);
  }

  {
    test::TestSupercell ts(primclex, Eigen::Vector3i(4, 3, 2).asDiagonal());
    local_pair_clusters_test(ts, phenom);
  }
}
