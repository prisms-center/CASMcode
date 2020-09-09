#include "gtest/gtest.h"

/// What is being tested:
#include "casm/clex/ConfigEnumRandomOccupations.hh"

/// What is being used to test it:
#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/app/enum.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/completer/Handlers.hh"
#include "casm/database/Database.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"

using namespace CASM;

TEST(ConfigEnumRandomOccupationsTest, Test1) {

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
    EXPECT_EQ(n_configs, std::distance(e.begin(), e.end()));
  }


  {
    jsonParser json;
    json["existing_only"] = false;
    json["max"] = 4;

    ScelEnum e(primclex, json);
    for(const auto &scel : e) {
      (void) scel;
    }
    EXPECT_EQ(primclex.generic_db<Supercell>().size(), 14);
  }

  {
    Completer::EnumOption enum_opt;
    jsonParser json;
    json["n_config"] = 200;
    ConfigEnumRandomOccupations::run(primclex, json, enum_opt, nullptr);
    EXPECT_GT(primclex.generic_db<Configuration>().size(), 150);
  }

}

TEST(ConfigEnumRandomOccupationsTest, ConfigEnumRandomOccupationsRunTest) {

  // create a project
  test::FCCTernaryProj proj;
  proj.check_init();

  // construct PrimClex
  PrimClex primclex(proj.dir, null_log());

  {
    Completer::EnumOption opt;
    parse_args(opt, "casm enum --method ScelEnum --max 4", primclex);
    ScelEnum::run(primclex, jsonParser(), opt, nullptr);
  }
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  primclex.generic_db<Supercell>().close();
  primclex.generic_db<Supercell>().open();
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);

  std::string cmd = R"(casm enum --method ConfigEnumRandomOccupations -a)";
  jsonParser kwargs;
  kwargs.put_obj();
  kwargs["n_config"] = 10;

  // --dry-run test
  {
    Completer::EnumOption opt;
    parse_args(opt, cmd + " --dry-run", primclex);
    ConfigEnumRandomOccupations::run(primclex, kwargs, opt, nullptr);
  }
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  EXPECT_GT(primclex.generic_db<Configuration>().size(), 0);
  primclex.generic_db<Configuration>().close();
  primclex.generic_db<Configuration>().open();
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  EXPECT_EQ(primclex.generic_db<Configuration>().size(), 0);

  {
    Completer::EnumOption opt;
    parse_args(opt, cmd, primclex);
    ConfigEnumRandomOccupations::run(primclex, kwargs, opt, nullptr);
  }
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  EXPECT_GT(primclex.generic_db<Configuration>().size(), 0);
  primclex.generic_db<Configuration>().close();
  primclex.generic_db<Configuration>().open();
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  EXPECT_GT(primclex.generic_db<Configuration>().size(), 0);

}
