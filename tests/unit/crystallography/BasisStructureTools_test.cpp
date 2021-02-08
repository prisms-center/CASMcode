#include "autotools.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

using namespace CASM;

TEST(BasicStructureToolsTest, MakeAsymmetricUnitTest) {
  EXPECT_EQ(xtal::make_asymmetric_unit(test::no_basis_prim()),
            std::set<std::set<Index>>());
  EXPECT_EQ(xtal::make_asymmetric_unit(test::FCC_ternary_prim()),
            std::set<std::set<Index>>({{0}}));
  EXPECT_EQ(xtal::make_asymmetric_unit(test::ZrO_prim()),
            std::set<std::set<Index>>({{0, 1}, {2, 3}}));
}
