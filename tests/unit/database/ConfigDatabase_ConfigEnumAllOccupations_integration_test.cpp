
#include "casm/app/ProjectBuilder.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/database/ConfigDatabaseTools_impl.hh"
#include "casm/database/Database.hh"
#include "casm/database/ScelDatabaseTools_impl.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

using namespace CASM;

TEST(ConfigDatabase_ConfigEnumAllOccupations_IntegrationTest, Test1) {
  // ConfigEnumAllOccupations: ZrO enumeration, w/ database

  auto shared_prim = std::make_shared<Structure const>(test::ZrO_prim());
  auto title = shared_prim->structure().title();
  auto project_settings = make_default_project_settings(*shared_prim, title);
  PrimClex primclex{project_settings, shared_prim};

  xtal::ScelEnumProps scel_enum_props{1, 5};
  ScelEnumByProps supercell_enumerator{shared_prim, scel_enum_props};
  for (auto const &supercell : supercell_enumerator) {
    supercell.set_primclex(&primclex);
    make_canonical_and_insert(supercell_enumerator, supercell,
                              primclex.db<Supercell>());
  }

  // order should not change as long as supercell comparison is not changed
  std::vector<Index> number_expected{3,  3,  4,  3,  16, 10, 10, 10, 10, 18,
                                     24, 24, 21, 15, 27, 27, 42, 24, 24, 21};

  std::vector<Index> number_enumerated;
  bool primitive_only = true;  //
  for (auto const &supercell : primclex.db<Supercell>()) {
    ConfigEnumAllOccupations enumerator{supercell};
    for (auto const &configuration : enumerator) {
      make_canonical_and_insert(enumerator, configuration,
                                primclex.db<Supercell>(),
                                primclex.db<Configuration>(), primitive_only);
    }
    auto range = primclex.db<Configuration>().scel_range(supercell.name());
    number_enumerated.push_back(std::distance(range.begin(), range.end()));
  }
  EXPECT_EQ(number_enumerated.size(), number_expected.size());
  EXPECT_EQ(number_enumerated, number_expected);
}
