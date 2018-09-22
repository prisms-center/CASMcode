#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/ConfigEnumRandomOccupations.hh"

/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/clex/ScelEnum.hh"
#include "casm/completer/Handlers.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(ConfigEnumRandomOccupationsTest)

BOOST_AUTO_TEST_CASE(Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.get_prim().lattice().vectors();

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
    for(const auto &scel : e) {
      (void) scel;
    }
    BOOST_CHECK_EQUAL(primclex.get_supercell_list().size(), 14);
  }

  {
    Completer::EnumOption enum_opt;
    jsonParser json;
    json["n_configs"] = 200;
    ConfigEnumRandomOccupations::run(primclex, json, enum_opt);
    for(const auto &scel : primclex.get_supercell_list()) {
      BOOST_CHECK_EQUAL(scel.get_config_list().size() > 0, true);
    }
  }

}

BOOST_AUTO_TEST_SUITE_END()
