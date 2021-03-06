#include "gtest/gtest.h"

/// What is being tested:
#include "casm/clex/ConfigEnumRandomOccupations.hh"

/// What is being used to test it:
#include "App/TestEnumeratorInterface.hh"
#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/app/enum.hh"
#include "casm/app/enum/methods/ConfigEnumRandomOccupationsInterface.hh"
#include "casm/app/enum/methods/ScelEnumInterface.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/completer/Handlers.hh"
#include "casm/database/Database.hh"
#include "casm/database/ScelDatabaseTools_impl.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"

using namespace CASM;

TEST(ConfigEnumRandomOccupationsTest, Test1) {
  ScopedNullLogging logging;

  test::FCCTernaryProj proj;
  proj.check_init();
  PrimClex primclex(proj.dir);
  auto &supercell_db = primclex.db<Supercell>();

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  Supercell scel =
      Supercell(&primclex, Lattice{2. * a, 2. * b, 2. * c}).canonical_form();

  Index n_configs = 100;
  MTRand mtrand;

  {
    ConfigEnumRandomOccupations e(scel, n_configs, mtrand);
    EXPECT_EQ(n_configs, std::distance(e.begin(), e.end()));
  }

  {
    xtal::ScelEnumProps enum_props(1, 5);
    ScelEnumByProps enumerator{primclex.shared_prim(), enum_props};
    for (const auto &supercell : enumerator) {
      make_canonical_and_insert(enumerator, supercell, supercell_db);
    }
    EXPECT_EQ(primclex.generic_db<Supercell>().size(), 14);
  }

  {
    jsonParser json_options;
    json_options["n_config"] = 200;
    std::string cli_str = "casm enum --method ConfigEnumRandomOccupations -a";
    test::run_enum_interface<ConfigEnumRandomOccupationsInterface>(
        cli_str, primclex, json_options);
    EXPECT_GT(primclex.generic_db<Configuration>().size(), 150);
  }
}
