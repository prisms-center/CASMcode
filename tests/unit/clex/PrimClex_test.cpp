#include "gtest/gtest.h"

/// What is being tested:
#include "casm/clex/PrimClex.hh"

/// What is being used to test it:

#include "casm/app/ProjectBuilder.hh"
#include "casm/crystallography/Structure.hh"
#include "Common.hh"
#include "FCCTernaryProj.hh"

using namespace CASM;

TEST(PrimClexTest, Basics) {

  auto shared_prim = std::make_shared<Structure const>(test::FCC_ternary_prim());
  auto const &prim = shared_prim->structure();
  EXPECT_EQ(prim.basis().size(), 1);

  // Construct PrimClex from project settings & prim
  auto project_settings = make_default_project_settings(prim, prim.title());
  PrimClex primclex {project_settings, shared_prim, null_log()};
  EXPECT_EQ(primclex.prim().basis().size(), 1);

}
