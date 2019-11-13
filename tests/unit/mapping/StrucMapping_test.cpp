#include "gtest/gtest.h"
#include "autotools.hh"

/// What is being tested:
#include "casm/crystallography/StrucMapping.hh"

/// What is being used to test it:
#include "crystallography/TestStructures.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/SimpleStrucMapCalculator.hh"

using namespace CASM;


/** PRIM1 *****************************
Face-centered Cubic (FCC, cF)
1.0
0 2.0 2.0
2.0 0 2.0
2.0 2.0 0
1
D
0.00 0.00 0.00 A B C
***************************************/

void sym_mapping_test(xtal::BasicStructure<xtal::Site> struc, Index N) {
  for(Index i = 0; i < struc.basis().size(); ++i)
    struc.set_occ(i, 0);

  xtal::SimpleStructure sstruc = xtal::to_simple_structure(struc);
  xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc)), 0.5, 0);

  auto sym_set = mapper.map_deformed_struc_impose_lattice(sstruc, xtal::Lattice(sstruc.lat_column_mat), xtal::StrucMapping::big_inf(), 1e-3);
  EXPECT_EQ(sym_set.size(), N);
}



TEST(SymMappingTest, FCCTernaryPrim) {
  // Read in test PRIM and run tests
  sym_mapping_test(test::FCC_ternary_prim(), 48);
}

TEST(SymMappingTest, ZrOPrim) {
  // Read in test PRIM and run tests
  sym_mapping_test(test::ZrO_prim(), 24);
}


