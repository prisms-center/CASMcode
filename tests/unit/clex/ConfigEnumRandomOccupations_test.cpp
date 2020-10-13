#include "gtest/gtest.h"

/// What is being tested:
#include "casm/clex/ConfigEnumRandomOccupations.hh"

/// What is being used to test it:
#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "clex/TestClexEnumerators.hh"
#include "casm/app/enum.hh"
#include "casm/app/enum/methods/ConfigEnumRandomOccupationsInterface.hh"
#include "casm/app/enum/methods/ScelEnumInterface.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/completer/Handlers.hh"
#include "casm/database/Database.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"

using namespace CASM;

TEST(ConfigEnumRandomOccupationsTest, Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  ScopedNullLogging logging;
  PrimClex primclex(proj.dir);

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  Supercell scel = Supercell(&primclex, Lattice {2.*a, 2.*b, 2.*c}).canonical_form();

  Index n_configs = 100;
  MTRand mtrand;

  {
    ConfigEnumRandomOccupations e(scel, n_configs, mtrand);
    EXPECT_EQ(n_configs, std::distance(e.begin(), e.end()));
  }


  {
    xtal::ScelEnumProps enum_props(1, 5);
    ScelEnumByProps e(primclex.shared_prim(), enum_props);
    for(const auto &scel : e) {
      (void) scel;
    }
    EXPECT_EQ(primclex.generic_db<Supercell>().size(), 14);
  }

  {
    jsonParser json_options;
    json_options["n_config"] = 200;
    std::string cli_str = "casm enum --method ConfigEnumRandomOccupations -a";
    test::run_enum_interface<ConfigEnumRandomOccupationsInterface>(cli_str, primclex, json_options);
    EXPECT_GT(primclex.generic_db<Configuration>().size(), 150);
  }

}
