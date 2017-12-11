#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clusterography/ClusterSpecsParser_impl.hh"

/// What is being used to test it:
#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/kinetics/DiffusionTransformation_impl.hh"
#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"
#include "casm/casm_io/jsonFile.hh"
#include "casm/completer/Handlers.hh"
#include "FCCTernaryProj.hh"
#include "ZrOProj.hh"
#include "TestOrbits.hh"
#include "TestSupercell.hh"
#include "TestConfiguration.hh"

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
               PrimPeriodicSymCompare<IntegralCluster>(primclex),
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
        jsonFile("tests/unit/kinetics/FCCTernary_bspecs_0.json"),
        2, 4) {
      // counts only check for consistency, not accuracy
      BOOST_CHECK_EQUAL(orbits.size(), 31);
      BOOST_CHECK_EQUAL(diff_trans_orbits.size(), 516);
    }
  };

  struct TestOrbits1 : test::TestPrimPeriodicDiffusionTransformationOrbits {
    TestOrbits1(const PrimClex &primclex) :
      test::TestPrimPeriodicDiffusionTransformationOrbits(
        primclex,
        jsonFile("tests/unit/kinetics/ZrO_bspecs_0.json"),
        2, 4) {
      // counts only check for consistency, not accuracy
      BOOST_CHECK_EQUAL(orbits.size(), 74);
      BOOST_CHECK_EQUAL(diff_trans_orbits.size(), 181);
    }
  };

  void local_point_clusters_test(
    const test::TestSupercell &ts,
    const Kinetics::DiffusionTransformation &phenom) {

    IntegralCluster phenom_sites = phenom.cluster();

    std::vector<PermuteIterator> dt_permute_group = phenom.invariant_subgroup(ts.scel);
    SymGroup dt_group = make_sym_group(dt_permute_group);

    OrbitGenerators<WithinScelIntegralClusterOrbit> generators(
      dt_group, ts.within_scel_sym_compare);

    for(Index l = 0; l < ts.scel.num_sites(); ++l) {
      IntegralCluster el(ts.scel.primclex().prim());
      el.elements().push_back(ts.scel.uccoord(l));
      generators.insert(el);
    }

    std::vector<WithinScelIntegralClusterOrbit> local_point_orbits;
    generators.make_orbits(std::back_inserter(local_point_orbits));

    Index count = 0;
    for(const auto &orbit : local_point_orbits) {
      count += orbit.size();
    }

    BOOST_CHECK_EQUAL(count, ts.scel.num_sites());

    FullOrbitPrinter<IntegralCluster> orbit_printer(6, '\n', CART);
    print_clust(local_point_orbits.begin(), local_point_orbits.end(), ts.scel.primclex().log(), orbit_printer);

  }

  void local_pair_clusters_test(
    const test::TestSupercell &ts,
    const Kinetics::DiffusionTransformation &phenom) {

    IntegralCluster phenom_sites = phenom.cluster();

    std::vector<PermuteIterator> dt_permute_group = phenom.invariant_subgroup(ts.scel);
    SymGroup dt_group = make_sym_group(dt_permute_group);

    OrbitGenerators<WithinScelIntegralClusterOrbit> generators(
      dt_group, ts.within_scel_sym_compare);

    for(Index l = 0; l < ts.scel.num_sites(); ++l) {
      for(Index m = 0; m < ts.scel.num_sites(); ++m) {
        if(l < m) {
          IntegralCluster el(ts.scel.primclex().prim());
          el.elements().push_back(ts.scel.uccoord(l));
          el.elements().push_back(ts.scel.uccoord(m));
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

    BOOST_CHECK_EQUAL(count, ts.scel.num_sites() * (ts.scel.num_sites() - 1) / 2);

    FullOrbitPrinter<IntegralCluster> orbit_printer(6, '\n', CART);
    print_clust(local_pair_orbits.begin(), local_pair_orbits.end(), ts.scel.primclex().log(), orbit_printer);

  }

}

BOOST_AUTO_TEST_SUITE(ClusterSpecsParserTest)

BOOST_AUTO_TEST_CASE(PrimPeriodicClustersByMaxLengthTest) {

  /// Make test project
  test::FCCTernaryProj proj;
  proj.check_init();
  PrimClex primclex(proj.dir, null_log());
  BOOST_CHECK_EQUAL(true, true);

  auto &log = null_log();
  {
    auto parser = PrimParser(primclex, R"({
      "orbit_branch_specs": {
        "2": {"max_length": 3.}
      }
    })");

    BOOST_CHECK_EQUAL(parser->valid(), true);
    BOOST_CHECK_EQUAL(parser->orbit_branch_specs().warning.size(), 0);

    BOOST_CHECK_EQUAL(parser->max_branch(), 2);
    BOOST_CHECK_EQUAL(almost_equal(parser->max_length(2), 3.0), true);
    BOOST_CHECK_THROW(parser->max_length(1), std::invalid_argument);
    BOOST_CHECK_THROW(parser->max_length(3), std::invalid_argument);
    BOOST_CHECK_EQUAL(parser->custom_generators().elements.size(), 0);

    log << parser->report() << std::endl;
  }

  {
    auto parser = PrimParser(primclex, R"({})");

    BOOST_CHECK_EQUAL(parser->valid(), true);
    BOOST_CHECK_EQUAL(parser->orbit_branch_specs().warning.size(), 0);

    BOOST_CHECK_EQUAL(parser->max_branch(), 1);
    BOOST_CHECK_EQUAL(parser->custom_generators().elements.size(), 0);

    log << parser->report() << std::endl;
  }

  {
    auto parser = PrimParser(primclex, R"({"orbit_branch_specs": []})");

    BOOST_CHECK_EQUAL(parser->valid(), false);
    BOOST_CHECK_EQUAL(parser->orbit_branch_specs().error.size(), 1);

    log << parser->report() << std::endl;
  }

  {
    auto parser = PrimParser(primclex, R"({
      "orbit_branch_specs": {
        "1": {"max_length": 3.},
        "2": {"max_length": 3.}
      }
    })");

    BOOST_CHECK_EQUAL(parser->valid(), true);
    BOOST_CHECK_EQUAL(parser->orbit_branch_specs().warning.size(), 1);

    log << parser->report() << std::endl;
  }

  {
    auto parser = PrimParser(primclex, R"({
      "orbit_branch_specs": {
        "3": {"max_length": 3.}
      }
    })");

    BOOST_CHECK_EQUAL(parser->valid(), false);
    BOOST_CHECK_EQUAL(parser->orbit_branch_specs().warning.size(), 0);
    BOOST_CHECK_EQUAL(parser->orbit_branch_specs().error.size(), 1);

    log << parser->report() << std::endl;
  }

  {
    auto parser = PrimParser(primclex, R"({
      "orbit_branch_specs": {
        "2": {"max_length": 3.},
        "3": {"max_length": 4.}
      }
    })");

    BOOST_CHECK_EQUAL(parser->valid(), false);
    BOOST_CHECK_EQUAL(parser->orbit_branch_specs().warning.size(), 0);
    BOOST_CHECK_EQUAL(parser->orbit_branch_specs().error.size(), 1);

    log << parser->report() << std::endl;
  }

  {
    auto parser = PrimParser(primclex, R"({
      "orbit_branch_specs": {
        "2": {"max_length": 4.},
        "3": {"max_length": 4.}
      }
    })");

    BOOST_CHECK_EQUAL(parser->valid(), true);
    BOOST_CHECK_EQUAL(parser->orbit_branch_specs().warning.size(), 0);
    BOOST_CHECK_EQUAL(parser->orbit_branch_specs().error.size(), 0);

    log << parser->report() << std::endl;
  }

  {
    auto parser = PrimParser(primclex, R"({
      "orbit_branch_specs": {
        "2": {"max_length": 4., "something_else": 1},
        "3": {"max_length": 4.},
        "_3": {"max_length": 5.}
      }
    })");

    BOOST_CHECK_EQUAL(parser->valid(), true);
    BOOST_CHECK_EQUAL(parser->orbit_branch_specs().warning.size(), 2);
    BOOST_CHECK_EQUAL(parser->orbit_branch_specs().error.size(), 0);

    log << parser->report() << std::endl;
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

    BOOST_CHECK_EQUAL(parser->valid(), true);
    BOOST_CHECK_EQUAL(parser->orbit_specs().warning.size(), 0);
    BOOST_CHECK_EQUAL(parser->orbit_specs().error.size(), 0);
    BOOST_CHECK_EQUAL(parser->custom_generators().elements.size(), 3);

    log << parser->report() << std::endl;
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

    BOOST_CHECK_EQUAL(parser->valid(), true);
    BOOST_CHECK_EQUAL(parser->orbit_specs().warning.size(), 0);
    BOOST_CHECK_EQUAL(parser->orbit_specs().error.size(), 0);
    BOOST_CHECK_EQUAL(parser->custom_generators().elements.size(), 1);

    log << parser->report() << std::endl;
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

    BOOST_CHECK_EQUAL(parser->valid(), true);
    BOOST_CHECK_EQUAL(parser->orbit_specs().warning.size(), 0);
    BOOST_CHECK_EQUAL(parser->orbit_specs().error.size(), 0);
    BOOST_CHECK_EQUAL(parser->custom_generators().elements.size(), 4);

    log << parser->report() << std::endl;
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

    BOOST_CHECK_EQUAL(parser->valid(), true);
    BOOST_CHECK_EQUAL(parser->orbit_specs().warning.size(), 1);
    BOOST_CHECK_EQUAL(parser->orbit_specs().error.size(), 0);
    BOOST_CHECK_EQUAL(parser->custom_generators().elements.size(), 4);

    log << parser->report() << std::endl;
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

    BOOST_CHECK_EQUAL(parser->valid(), false);
    BOOST_CHECK_EQUAL(parser->orbit_specs().error.size(), 1);

    log << parser->report() << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(LocalClustersByMaxLengthTest) {
  /// Make test project
  test::FCCTernaryProj proj;
  proj.check_init();
  Log &log = default_log();
  PrimClex primclex(proj.dir, null_log());
  BOOST_CHECK_EQUAL(true, true);

  TestOrbits0 to(primclex);

  jsonFile diff_trans_json {"tests/unit/kinetics/FCCTernary_diff_trans_0.json"};
  Completer::EnumOption enum_opt;
  enum_opt.desc();
  int success = Kinetics::DiffusionTransformationEnum::run(primclex, diff_trans_json, enum_opt);
  const auto &diff_trans_db = primclex.generic_db<Kinetics::PrimPeriodicDiffTransOrbit>();
  BOOST_CHECK_EQUAL(diff_trans_db.size(), 28);
  BOOST_CHECK_EQUAL(success, 0);

  PrototypePrinter<Kinetics::DiffusionTransformation> dt_orbit_printer(6, '\n', CART);

  {
    Index i = 0;
    for(const auto &orbit : diff_trans_db) {
      log << "\n --- " << orbit.name() << " ---" << std::endl;
      dt_orbit_printer(orbit, log, i++, diff_trans_db.size());
    }
  }



  // generate local orbits (w/ scel symmetry) using a 3x3x3 standard FCC supercell (so 3x3x3x4 sites)
  test::TestStandardFCCSupercell ts(primclex, 3, 3, 3);

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
    auto parser = LocalParser(ts.scel, R"({
      "standard": {
        "orbit_branch_specs": {
          "1": {"cutoff_radius": 10.}
        }
      }
    })");

    BOOST_CHECK_EQUAL(parser->valid(), true);
    BOOST_CHECK_EQUAL(parser->standard->max_branch, 1);
    BOOST_CHECK_EQUAL(parser->standard->max_length_including_phenomenal(), false);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->cutoff_radius(1), 10.), true);
    BOOST_CHECK_THROW(parser->standard->max_length(1), std::runtime_error);

    BOOST_CHECK_EQUAL(parser->all_warning().size(), 0);
    BOOST_CHECK_EQUAL(parser->all_error().size(), 0);

    log << parser->report() << std::endl;
  }

  // --- typical ---
  {
    auto parser = LocalParser(ts.scel, R"({
      "standard": {
        "orbit_branch_specs": {
          "1": {"cutoff_radius": 3.},
          "2": {"max_length": 12.0, "cutoff_radius": 3.},
          "3": {"max_length": 10.0, "cutoff_radius": 3.}
        }
      }
    })");

    BOOST_CHECK_EQUAL(parser->valid(), true);
    BOOST_CHECK_EQUAL(parser->standard->max_branch, 3);
    BOOST_CHECK_EQUAL(parser->standard->max_length_including_phenomenal(), false);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->cutoff_radius(1), 3.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->cutoff_radius(2), 3.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->cutoff_radius(3), 3.), true);
    BOOST_CHECK_THROW(parser->standard->max_length(1), std::runtime_error);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->max_length(2), 12.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->max_length(3), 10.), true);

    BOOST_CHECK_EQUAL(parser->all_warning().size(), 0);
    BOOST_CHECK_EQUAL(parser->all_error().size(), 0);

    log << parser->report() << std::endl;
  }

  // --- typical, warning for max_length on brach 1 w/out max_length_including_phenomenal option ---
  {
    auto parser = LocalParser(ts.scel, R"({
      "standard": {
        "orbit_branch_specs": {
          "1": {"max_length": 12.0, "cutoff_radius": 3.},
          "2": {"max_length": 12.0, "cutoff_radius": 3.},
          "3": {"max_length": 10.0, "cutoff_radius": 3.}
        }
      }
    })");

    BOOST_CHECK_EQUAL(parser->valid(), true);
    BOOST_CHECK_EQUAL(parser->standard->warning.size(), 1);
    BOOST_CHECK_EQUAL(parser->standard->error.size(), 0);
    BOOST_CHECK_EQUAL(parser->all_warning().size(), 1);
    BOOST_CHECK_EQUAL(parser->all_error().size(), 0);

    log << parser->report() << std::endl;
  }

  // --- typical, error for nonincreasing ---
  {
    auto parser = LocalParser(ts.scel, R"({
      "standard": {
        "orbit_branch_specs": {
          "1": {"cutoff_radius": 3.},
          "2": {"max_length": 12.0, "cutoff_radius": 4.},
          "3": {"max_length": 13.0, "cutoff_radius": 3.}
        }
      }
    })");

    BOOST_CHECK_EQUAL(parser->valid(), false);
    BOOST_CHECK_EQUAL(parser->standard->warning.size(), 0);
    BOOST_CHECK_EQUAL(parser->standard->error.size(), 2);
    BOOST_CHECK_EQUAL(parser->all_warning().size(), 0);
    BOOST_CHECK_EQUAL(parser->all_error().size(), 2);

    log << parser->report() << std::endl;
  }

  // --- typical, w/ max_length_including_phenomenal ---
  {
    auto parser = LocalParser(ts.scel, R"({
      "standard": {
        "orbit_branch_specs": {
          "max_length_including_phenomenal": true,
          "1": {"max_length": 12.0, "cutoff_radius": 3.},
          "2": {"max_length": 12.0, "cutoff_radius": 3.},
          "3": {"max_length": 10.0, "cutoff_radius": 3.}
        }
      }
    })");

    BOOST_CHECK_EQUAL(parser->valid(), true);
    BOOST_CHECK_EQUAL(parser->standard->max_branch, 3);
    BOOST_CHECK_EQUAL(parser->standard->max_length_including_phenomenal(), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->cutoff_radius(1), 3.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->cutoff_radius(2), 3.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->cutoff_radius(3), 3.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->max_length(1), 12.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->max_length(2), 12.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->max_length(3), 10.), true);

    BOOST_CHECK_EQUAL(parser->all_warning().size(), 0);
    BOOST_CHECK_EQUAL(parser->all_error().size(), 0);

    log << parser->report() << std::endl;
  }

  // --- typical, w/ max_length_including_phenomenal && max_length="inf" (all clusters of branch in neighborhood) ---
  {
    auto parser = LocalParser(ts.scel, R"({
      "standard": {
        "orbit_branch_specs": {
          "max_length_including_phenomenal": true,
          "1": {"max_length": "inf", "cutoff_radius": 3.},
          "2": {"max_length": "inf", "cutoff_radius": 3.},
          "3": {"max_length": "inf", "cutoff_radius": 3.}
        }
      }
    })");

    BOOST_CHECK_EQUAL(parser->valid(), true);
    BOOST_CHECK_EQUAL(parser->standard->max_branch, 3);
    BOOST_CHECK_EQUAL(parser->standard->max_length_including_phenomenal(), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->cutoff_radius(1), 3.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->cutoff_radius(2), 3.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->cutoff_radius(3), 3.), true);
    BOOST_CHECK_EQUAL(parser->standard->max_length(1), std::numeric_limits<double>::infinity());
    BOOST_CHECK_EQUAL((3.0 < parser->standard->max_length(1)), true);
    BOOST_CHECK_EQUAL(parser->standard->max_length(2), std::numeric_limits<double>::infinity());
    BOOST_CHECK_EQUAL(parser->standard->max_length(3), std::numeric_limits<double>::infinity());

    BOOST_CHECK_EQUAL(parser->all_warning().size(), 0);
    BOOST_CHECK_EQUAL(parser->all_error().size(), 0);

    log << parser->report() << std::endl;
  }

  // --- w/ max_length_including_phenomenal, error for nonincreasing ---
  {
    auto parser = LocalParser(ts.scel, R"({
      "standard": {
        "orbit_branch_specs": {
          "max_length_including_phenomenal": true,
          "1": {"max_length": 13.0, "cutoff_radius": 3.},
          "2": {"max_length": 12.0, "cutoff_radius": 4.},
          "3": {"max_length": 13.0, "cutoff_radius": 3.}
        }
      }
    })");

    BOOST_CHECK_EQUAL(parser->valid(), false);
    BOOST_CHECK_EQUAL(parser->standard->warning.size(), 0);
    BOOST_CHECK_EQUAL(parser->standard->error.size(), 2);
    BOOST_CHECK_EQUAL(parser->all_warning().size(), 0);
    BOOST_CHECK_EQUAL(parser->all_error().size(), 2);

    log << parser->report() << std::endl;
  }

  // --- custom orbit_branch_specs, equivalence_type=prim (default)  ---
  {
    auto parser = LocalParser(ts.scel, R"({
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
            "specie_traj" : [
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

    BOOST_CHECK_EQUAL(parser->valid(), true);
    BOOST_CHECK_EQUAL(parser->standard->max_branch, 2);
    BOOST_CHECK_EQUAL(parser->standard->max_length_including_phenomenal(), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->cutoff_radius(1), 3.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->cutoff_radius(2), 3.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->max_length(1), 4.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->max_length(2), 4.), true);

    BOOST_CHECK_EQUAL(parser->custom->data.size(), 1);
    BOOST_CHECK_EQUAL(parser->custom->data[0].orbit_branch_specs->max_branch, 3);
    BOOST_CHECK_EQUAL(parser->custom->data[0].orbit_branch_specs->max_length_including_phenomenal(), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->custom->data[0].orbit_branch_specs->cutoff_radius(1), 6.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->custom->data[0].orbit_branch_specs->cutoff_radius(2), 6.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->custom->data[0].orbit_branch_specs->cutoff_radius(3), 6.), true);
    BOOST_CHECK_EQUAL(parser->custom->data[0].orbit_branch_specs->max_length(1), std::numeric_limits<double>::infinity());
    BOOST_CHECK_EQUAL(parser->custom->data[0].orbit_branch_specs->max_length(2), std::numeric_limits<double>::infinity());
    BOOST_CHECK_EQUAL(parser->custom->data[0].orbit_branch_specs->max_length(3), std::numeric_limits<double>::infinity());

    // checks comparing input phenomenal cluster to test cluster
    // the input custom phenomenal cluster is to.diff_trans_orbits[0].prototype(), so
    // only clusters in orbit to.diff_trans_orbits[0] should be found
    {
      Index linear_orbit_index = 0;
      for(const auto &orbit : to.diff_trans_orbits) {
        for(const auto &equiv : orbit) {
          auto find_res = parser->find(equiv);
          auto it = find_res.first;
          auto op = find_res.second;
          BOOST_CHECK_EQUAL((it != parser->custom->data.end()), (linear_orbit_index == 0));

          // check op*phenom == equiv
          if(it != parser->custom->data.end()) {
            auto &f = *it->phenom->prim_sym_compare;
            BOOST_CHECK_EQUAL(f.equal(f.prepare(copy_apply(op, *it->phenom->phenom)), f.prepare(equiv)), true);
          }

          int expected_max_branch = (linear_orbit_index == 0 ? 3 : 2);
          BOOST_CHECK_EQUAL(parser->max_branch(it), expected_max_branch);
          BOOST_CHECK_EQUAL(parser->max_length_including_phenomenal(it), true);

          for(int branch = 1; branch <= parser->max_branch(it); ++branch) {
            double expected_cutoff_radius = (linear_orbit_index == 0 ? 6. : 3.);
            BOOST_CHECK_EQUAL(parser->cutoff_radius(it, branch), expected_cutoff_radius);

            double expected_max_length = (linear_orbit_index == 0 ? std::numeric_limits<double>::infinity() : 4.);
            BOOST_CHECK_EQUAL(parser->max_length(it, branch), expected_max_length);
          }
        }
        linear_orbit_index++;
      }
    }

    BOOST_CHECK_EQUAL(parser->all_warning().size(), 0);
    BOOST_CHECK_EQUAL(parser->all_error().size(), 0);

    log << parser->report() << std::endl;
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
    //    BOOST_CHECK_EQUAL(parser->valid(), true);
    //    BOOST_CHECK_EQUAL(parser->standard->max_branch, 2);
    //    BOOST_CHECK_EQUAL(parser->standard->max_length_including_phenomenal(), true);
    //    BOOST_CHECK_EQUAL(almost_equal(parser->standard->cutoff_radius(1), 3.), true);
    //    BOOST_CHECK_EQUAL(almost_equal(parser->standard->cutoff_radius(2), 3.), true);
    //    BOOST_CHECK_EQUAL(almost_equal(parser->standard->max_length(1), 4.), true);
    //    BOOST_CHECK_EQUAL(almost_equal(parser->standard->max_length(2), 4.), true);
    //
    //    BOOST_CHECK_EQUAL(parser->custom->data.size(), 1);
    //    BOOST_CHECK_EQUAL(parser->custom->data[0].orbit_branch_specs->max_branch, 3);
    //    BOOST_CHECK_EQUAL(parser->custom->data[0].orbit_branch_specs->max_length_including_phenomenal(), true);
    //    BOOST_CHECK_EQUAL(almost_equal(parser->custom->data[0].orbit_branch_specs->cutoff_radius(1), 6.), true);
    //    BOOST_CHECK_EQUAL(almost_equal(parser->custom->data[0].orbit_branch_specs->cutoff_radius(2), 6.), true);
    //    BOOST_CHECK_EQUAL(almost_equal(parser->custom->data[0].orbit_branch_specs->cutoff_radius(3), 6.), true);
    //    BOOST_CHECK_EQUAL(parser->custom->data[0].orbit_branch_specs->max_length(1), std::numeric_limits<double>::infinity());
    //    BOOST_CHECK_EQUAL(parser->custom->data[0].orbit_branch_specs->max_length(2), std::numeric_limits<double>::infinity());
    //    BOOST_CHECK_EQUAL(parser->custom->data[0].orbit_branch_specs->max_length(3), std::numeric_limits<double>::infinity());
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
    //          BOOST_CHECK_EQUAL((it != parser->custom->data.end()), (linear_orbit_index == 0));
    //
    //          // check op*phenom == equiv
    //          if(it != parser->custom->data.end()) {
    //            auto &f = *it->phenom->prim_sym_compare;
    //            BOOST_CHECK_EQUAL(f.equal(f.prepare(copy_apply(op, *it->phenom->phenom)), f.prepare(equiv)), true);
    //          }
    //
    //          int expected_max_branch = (linear_orbit_index == 0 ? 3 : 2);
    //          BOOST_CHECK_EQUAL(parser->max_branch(it), expected_max_branch);
    //          BOOST_CHECK_EQUAL(parser->max_length_including_phenomenal(it), true);
    //
    //          for(int branch = 1; branch <= parser->max_branch(it); ++branch) {
    //            double expected_cutoff_radius = (linear_orbit_index == 0 ? 6. : 3.);
    //            BOOST_CHECK_EQUAL(parser->cutoff_radius(it, branch), expected_cutoff_radius);
    //
    //            double expected_max_length = (linear_orbit_index == 0 ? std::numeric_limits<double>::infinity() : 4.);
    //            BOOST_CHECK_EQUAL(parser->max_length(it, branch), expected_max_length);
    //          }
    //        }
    //        linear_orbit_index++;
    //      }
    //    }
    //
    //    BOOST_CHECK_EQUAL(parser->all_warning().size(), 0);
    //    BOOST_CHECK_EQUAL(parser->all_error().size(), 0);
    //
    //    log << parser->report() << std::endl;
  }

  // --- custom orbit_specs, equivalence_type=prim (default)  ---
  {
    auto parser = LocalParser(ts.scel, R"({
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
            "specie_traj" : [
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

    BOOST_CHECK_EQUAL(parser->valid(), true);
    BOOST_CHECK_EQUAL(parser->standard->max_branch, 2);
    BOOST_CHECK_EQUAL(parser->standard->max_length_including_phenomenal(), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->cutoff_radius(1), 3.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->cutoff_radius(2), 3.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->max_length(1), 4.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->max_length(2), 4.), true);

    BOOST_CHECK_EQUAL(parser->custom->data.size(), 1);
    BOOST_CHECK_EQUAL(parser->custom->data[0].orbit_specs->prototypes.size(), 1);
    BOOST_CHECK_EQUAL(parser->custom->data[0].orbit_specs->prototypes[0].include_subclusters, 1);
    BOOST_CHECK_EQUAL(parser->custom->data[0].orbit_specs->prototypes[0].cluster.size(), 1);

    ScelPeriodicDiffTransSymCompare dt_scel_sym_compare(ts.scel);
    // checks comparing input phenomenal cluster to test cluster
    // the input custom phenomenal cluster is to.diff_trans_orbits[0].prototype(), so
    // only clusters in orbit to.diff_trans_orbits[0] should be found
    {
      Index linear_orbit_index = 0;
      for(const auto &orbit : to.diff_trans_orbits) {
        Index equiv_index = 0;
        for(const auto &equiv : orbit) {

          auto find_res = parser->find(equiv);
          auto it = find_res.first;
          auto op = find_res.second;

          // only linear_orbit_index == 0 should match custom orbit specs
          BOOST_CHECK_EQUAL((it != parser->custom->data.end()), (linear_orbit_index == 0));

          // generate local orbit
          if(it != parser->custom->data.end()) {

            // check op*phenom == equiv
            auto &f = *it->phenom->prim_sym_compare;
            BOOST_CHECK_EQUAL(f.equal(f.prepare(copy_apply(op, *it->phenom->phenom)), f.prepare(equiv)), true);

            std::vector<PermuteIterator> dt_permute_group = equiv.invariant_subgroup(ts.scel);
            SymGroup dt_group = make_sym_group(dt_permute_group);
            OrbitGenerators<ScelPeriodicIntegralClusterOrbit> test_gen(dt_group, ts.scel_sym_compare);
            parser->insert_custom_generators(find_res, test_gen);

            if(linear_orbit_index == 0) {

              BOOST_CHECK_EQUAL(dt_group.size(), 8);

              // 2 elements (because default==include_subclusters -> null cluster)
              BOOST_CHECK_EQUAL(test_gen.elements.size(), 2);

              auto gen_it = test_gen.elements.begin();
              {
                // orbit of null cluster -> 1 equiv
                ScelPeriodicIntegralClusterOrbit orbit(*gen_it++, dt_group, ts.scel_sym_compare);
                BOOST_CHECK_EQUAL(orbit.size(), 1);
              }

              {
                // should generate 4 equivalents (there are 4NN of both sites in 1NN pair cluster)
                ScelPeriodicIntegralClusterOrbit orbit(*gen_it++, dt_group, ts.scel_sym_compare);
                BOOST_CHECK_EQUAL(orbit.size(), 4);
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

    BOOST_CHECK_EQUAL(parser->all_warning().size(), 0);
    BOOST_CHECK_EQUAL(parser->all_error().size(), 0);

    log << parser->report() << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(LocalClustersByMaxLengthTest_Tet) {
  /// Make test project
  test::FCCTernaryProj proj;
  proj.check_init();
  Log &log = default_log();
  PrimClex primclex(proj.dir, null_log());
  BOOST_CHECK_EQUAL(true, true);

  TestOrbits0 to(primclex);


  // generate local orbits (w/ scel symmetry) using a 4x3x3 standard FCC supercell (so 4x3x3x4 sites)
  test::TestStandardFCCSupercell ts(primclex, 4, 3, 3);

  // --- custom orbit_specs, equivalence_type=scel, non-cubic supercell  ---
  {
    auto parser = LocalParser(ts.scel, R"({
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
            "specie_traj" : [
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
            "specie_traj" : [
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

    log << parser->report() << std::endl;

    BOOST_CHECK_EQUAL(parser->valid(), true);
    BOOST_CHECK_EQUAL(parser->standard->max_branch, 2);
    BOOST_CHECK_EQUAL(parser->standard->max_length_including_phenomenal(), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->cutoff_radius(1), 3.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->cutoff_radius(2), 3.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->max_length(1), 4.), true);
    BOOST_CHECK_EQUAL(almost_equal(parser->standard->max_length(2), 4.), true);

    BOOST_CHECK_EQUAL(parser->custom->data.size(), 2);
    BOOST_CHECK_EQUAL(parser->custom->data[0].orbit_specs->prototypes.size(), 1);
    BOOST_CHECK_EQUAL(parser->custom->data[0].orbit_specs->prototypes[0].include_subclusters, 1);
    BOOST_CHECK_EQUAL(parser->custom->data[0].orbit_specs->prototypes[0].cluster.size(), 1);
    BOOST_CHECK_EQUAL(parser->custom->data[1].orbit_specs->prototypes.size(), 1);
    BOOST_CHECK_EQUAL(parser->custom->data[1].orbit_specs->prototypes[0].include_subclusters, 1);
    BOOST_CHECK_EQUAL(parser->custom->data[1].orbit_specs->prototypes[0].cluster.size(), 1);

    ScelPeriodicDiffTransSymCompare dt_scel_sym_compare(ts.scel);
    // checks comparing input phenomenal cluster to test cluster
    // the input custom phenomenal cluster is to.diff_trans_orbits[0].prototype(), so
    // 4/6 clusters in orbit to.diff_trans_orbits[0] should match the first custom phenom
    // 2/6 should match the other custom phenom

    Printer<IntegralCluster> cluster_printer(6, '\n', CART);

    //    log << "check # matches" << std::endl;
    {
      std::vector<Index> count = {0, 0};
      BOOST_CHECK_EQUAL(to.diff_trans_orbits[0].size(), 6);
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
          BOOST_CHECK_EQUAL(f.equal(f.prepare(copy_apply(op, *it->phenom->phenom)), f.prepare(equiv)), true);

          std::vector<PermuteIterator> dt_permute_group = equiv.invariant_subgroup(ts.scel);
          SymGroup dt_group = make_sym_group(dt_permute_group);

          if(&*it == &parser->custom->data[0]) {
            count[0]++;
            BOOST_CHECK_EQUAL(dt_group.size(), 4);
          }
          else if(&*it == &parser->custom->data[1]) {
            count[1]++;
            BOOST_CHECK_EQUAL(dt_group.size(), 8);
          }

          //          log << "--- dt_group ---" << std::endl;
          //          Index op_index = 0;
          //          for(const auto& op: dt_group) {
          //            log << op_index++ << ":\n" << std::endl;
          //            log << "op.matrix(): \n" << op.matrix() << std::endl;
          //            log << "op.tau(): " << op.tau().transpose() << std::endl << std::endl;
          //          }

          OrbitGenerators<ScelPeriodicIntegralClusterOrbit> test_gen(dt_group, ts.scel_sym_compare);
          parser->insert_custom_generators(find_res, test_gen);

          // 2 elements (because default==include_subclusters -> null cluster)
          BOOST_CHECK_EQUAL(test_gen.elements.size(), 2);

          auto gen_it = test_gen.elements.begin();
          {
            // orbit of null cluster -> 1 equiv
            ScelPeriodicIntegralClusterOrbit orbit(*gen_it++, dt_group, ts.scel_sym_compare);
            BOOST_CHECK_EQUAL(orbit.size(), 1);
          }

          {
            // should generate 4 equivalents (the 4 NN of both sites in 1NN pair cluster)
            ScelPeriodicIntegralClusterOrbit orbit(*gen_it++, dt_group, ts.scel_sym_compare);
            BOOST_CHECK_EQUAL(orbit.size(), 4);
          }
        }
      }
      BOOST_CHECK_EQUAL(count[0], 4); // 4/6 equiv match first phenom
      BOOST_CHECK_EQUAL(count[1], 2); // 2/6 equiv match second phenom
    }

    // other orbits should not match (this is slow)
    log << "check # not matches" << std::endl;
    {
      for(Index i = 1; i < to.diff_trans_orbits.size(); ++i) {
        for(const auto &equiv : to.diff_trans_orbits[i]) {
          auto find_res = parser->find(equiv);
          auto it = find_res.first;
          auto op = find_res.second;
          BOOST_CHECK_EQUAL((it == parser->custom->data.end()), true);
        }
      }
    }

    BOOST_CHECK_EQUAL(parser->all_warning().size(), 0);
    BOOST_CHECK_EQUAL(parser->all_error().size(), 0);

    log << parser->report() << std::endl;
  }

}

BOOST_AUTO_TEST_CASE(LocalPointClustersByMaxLengthTest_FCCTernary) {
  /// Make test project
  test::FCCTernaryProj proj;
  proj.check_init();
  Log &log = default_log();
  PrimClex primclex(proj.dir, null_log());
  BOOST_CHECK_EQUAL(true, true);

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

BOOST_AUTO_TEST_CASE(LocalPointClustersByMaxLengthTest_ZrO) {
  /// Make test project
  test::ZrOProj proj;
  proj.check_init();
  Log &log = default_log();
  PrimClex primclex(proj.dir, null_log());
  BOOST_CHECK_EQUAL(true, true);

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

BOOST_AUTO_TEST_SUITE_END()
