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
#include "FCCTernaryProj.hh"

using namespace CASM;

std::unique_ptr<PrimPeriodicClustersByMaxLength>
prim_parser(const PrimClex &primclex, std::string specs) {
  return notstd::make_unique<PrimPeriodicClustersByMaxLength>(
           primclex,
           primclex.prim().factor_group(),
           PrimPeriodicSymCompare<IntegralCluster>(primclex),
           jsonParser::parse(specs),
           fs::path(),
           true);
}

BOOST_AUTO_TEST_SUITE(ClusterSpecsParserTest)

BOOST_AUTO_TEST_CASE(ClusterographyTest) {

  /// Make test project
  test::FCCTernaryProj proj;
  proj.check_init();
  PrimClex primclex(proj.dir, null_log());
  BOOST_CHECK_EQUAL(true, true);

  auto &log = null_log();
  {
    auto parser = prim_parser(primclex, R"({
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
    auto parser = prim_parser(primclex, R"({})");

    BOOST_CHECK_EQUAL(parser->valid(), true);
    BOOST_CHECK_EQUAL(parser->orbit_branch_specs().warning.size(), 0);

    BOOST_CHECK_EQUAL(parser->max_branch(), 1);
    BOOST_CHECK_EQUAL(parser->custom_generators().elements.size(), 0);

    log << parser->report() << std::endl;
  }

  {
    auto parser = prim_parser(primclex, R"({"orbit_branch_specs": []})");

    BOOST_CHECK_EQUAL(parser->valid(), false);
    BOOST_CHECK_EQUAL(parser->orbit_branch_specs().error.size(), 1);

    log << parser->report() << std::endl;
  }

  {
    auto parser = prim_parser(primclex, R"({
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
    auto parser = prim_parser(primclex, R"({
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
    auto parser = prim_parser(primclex, R"({
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
    auto parser = prim_parser(primclex, R"({
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
    auto parser = prim_parser(primclex, R"({
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
    auto parser = prim_parser(primclex, R"({
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
    auto parser = prim_parser(primclex, R"({
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
    auto parser = prim_parser(primclex, R"({
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
    auto parser = prim_parser(primclex, R"({
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
    auto parser = prim_parser(primclex, R"({
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

BOOST_AUTO_TEST_SUITE_END()
