#include "gtest/gtest.h"
#include "autotools.hh"

/// What is being tested:
#include "casm/crystallography/StrucMapping.hh"

/// What is being used to test it:
#include "crystallography/TestStructures.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/Adapter.hh"
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

  Eigen::Matrix3i T;
  T.setIdentity();
  T *= 2;

  xtal::SimpleStructure sstruc2 = to_superstructure(T, sstruc);

  // Check that we find 8 perfect mapping for a Vol8 non-primitive structure
  {
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc2)));
    xtal::LatticeNode tnode((xtal::Lattice(sstruc2.lat_column_mat)),
                            (xtal::Lattice(sstruc2.lat_column_mat)),
                            (xtal::Lattice(sstruc2.lat_column_mat)),
                            (xtal::Lattice(sstruc2.lat_column_mat)),
                            sstruc2.atom_info.size());

    auto trans_set = mapper.map_deformed_struc_impose_lattice_node(sstruc2, tnode, 0, xtal::StrucMapping::big_inf(), 1e-3);
    EXPECT_EQ(trans_set.size(), 8);
  }

  // Check for perfect mappings using the best-1 calling convention, without symmetry
  {
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(sstruc, xtal::Lattice(sstruc.lat_column_mat), 1, xtal::StrucMapping::big_inf(), -1e-3);
    EXPECT_EQ(sym_set.size(), N);
  }

  // Check for perfect mappings using the best-0 calling convention, without symmetry and with a positive min_cost
  // Store result as factor group of structure
  xtal::SymOpVector fgroup;
  {
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(sstruc, xtal::Lattice(sstruc.lat_column_mat), 0, xtal::StrucMapping::big_inf(), 1e-3);
    //pass
    EXPECT_EQ(sym_set.size(), N);
    fgroup = adapter::Adapter<xtal::SymOpVector, decltype(sym_set)>()(sym_set);
  }

  // Check for perfect mappings of primitive structure onto itself, using symmetry reduction of factor group from previous step.
  {
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc, fgroup)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(sstruc, xtal::Lattice(sstruc.lat_column_mat), 0, xtal::StrucMapping::big_inf(), 1e-3);
    //pass
    EXPECT_EQ(sym_set.size(), 1);
  }

  // Check for perfect mappings of non-primitive structure onto primitive, using symmetry reduction of factor group from previous step.
  {
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc, fgroup)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(sstruc2, xtal::Lattice(sstruc.lat_column_mat), 0, xtal::StrucMapping::big_inf(), 1e-3);
    EXPECT_EQ(sym_set.size(), 1);
  }

  // Check for perfect mappings of vol-8 non-primitive structure onto itself, using symmetry reduction of factor group from previous step.
  {
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc, fgroup)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(sstruc2, xtal::Lattice(sstruc.lat_column_mat), 0, xtal::StrucMapping::big_inf(), 1e-3);
    EXPECT_EQ(sym_set.size(), 8);
  }



}



TEST(SymMappingTest, FCCTernaryPrim) {
  // Read in test PRIM and run tests
  sym_mapping_test(test::FCC_ternary_prim(), 48);
}

TEST(SymMappingTest, ZrOPrim) {
  // Read in test PRIM and run tests
  sym_mapping_test(test::ZrO_prim(), 24);
}


