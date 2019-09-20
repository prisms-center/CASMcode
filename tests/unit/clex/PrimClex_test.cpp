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

  Structure prim(test::FCC_ternary_prim());
  EXPECT_EQ(prim.basis().size(), 1);


  // Construct from prim
  PrimClex primclex(prim, null_log());
  EXPECT_EQ(primclex.prim().basis().size(), 1);

}
