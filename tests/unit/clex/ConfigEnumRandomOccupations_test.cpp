#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/ConfigEnumRandomOccupations.hh"

/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/app/enum.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/database/Database.hh"
#include "casm/completer/Handlers.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(ConfigEnumRandomOccupationsTest)

BOOST_AUTO_TEST_CASE(Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  Supercell scel = Supercell(&primclex, Lattice {2.*a, 2.*b, 2.*c}).canonical_form();

  Index n_configs = 100;
  MTRand mtrand;

  {
    ConfigEnumRandomOccupations e(scel, n_configs, mtrand);
    BOOST_CHECK_EQUAL(n_configs, std::distance(e.begin(), e.end()));
  }


  {
    jsonParser json;
    json["existing_only"] = false;
    json["max"] = 4;

    ScelEnum e(primclex, json);
    for(const auto &scel : e) {}
    BOOST_CHECK_EQUAL(primclex.generic_db<Supercell>().size(), 14);
  }

  {
    Completer::EnumOption enum_opt;
    jsonParser json;
    json["n_config"] = 200;
    ConfigEnumRandomOccupations::run(primclex, json, enum_opt);
    BOOST_CHECK_EQUAL(primclex.generic_db<Configuration>().size() > 150, true);
  }

}

BOOST_AUTO_TEST_CASE(ConfigEnumRandomOccupationsRunTest) {

  // create a project
  test::FCCTernaryProj proj;
  proj.check_init();

  // construct PrimClex
  PrimClex primclex(proj.dir, null_log());

  {
    Completer::EnumOption opt;
    parse_args(opt, "casm enum --method ScelEnum --max 4", primclex);
    ScelEnum::run(primclex, jsonParser(), opt);
  }
  BOOST_CHECK_EQUAL(primclex.generic_db<Supercell>().size(), 13);
  primclex.generic_db<Supercell>().close();
  primclex.generic_db<Supercell>().open();
  BOOST_CHECK_EQUAL(primclex.generic_db<Supercell>().size(), 13);

  std::string cmd = R"(casm enum --method ConfigEnumRandomOccupations -a)";
  jsonParser kwargs;
  kwargs.put_obj();
  kwargs["n_config"] = 10;

  // --dry-run test
  {
    Completer::EnumOption opt;
    parse_args(opt, cmd + " --dry-run", primclex);
    ConfigEnumRandomOccupations::run(primclex, kwargs, opt);
  }
  BOOST_CHECK_EQUAL(primclex.generic_db<Supercell>().size(), 13);
  BOOST_CHECK_EQUAL(primclex.generic_db<Configuration>().size() > 0, true);
  primclex.generic_db<Configuration>().close();
  primclex.generic_db<Configuration>().open();
  BOOST_CHECK_EQUAL(primclex.generic_db<Supercell>().size(), 13);
  BOOST_CHECK_EQUAL(primclex.generic_db<Configuration>().size(), 0);

  {
    Completer::EnumOption opt;
    parse_args(opt, cmd, primclex);
    ConfigEnumRandomOccupations::run(primclex, kwargs, opt);
  }
  BOOST_CHECK_EQUAL(primclex.generic_db<Supercell>().size(), 13);
  BOOST_CHECK_EQUAL(primclex.generic_db<Configuration>().size() > 0, true);
  primclex.generic_db<Configuration>().close();
  primclex.generic_db<Configuration>().open();
  BOOST_CHECK_EQUAL(primclex.generic_db<Supercell>().size(), 13);
  BOOST_CHECK_EQUAL(primclex.generic_db<Configuration>().size() > 0, true);

}

BOOST_AUTO_TEST_SUITE_END()
