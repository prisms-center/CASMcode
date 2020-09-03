#include "gtest/gtest.h"

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/UnitCellCoord.hh"

#include "crystallography/TestStructures.hh" // for test::ZrO_prim

TEST(ExampleClusterographyUnitCellCoord, UnitCellCoordBasics) {

  using CASM::xtal::UnitCell;
  using CASM::xtal::UnitCellCoord;
  using CASM::xtal::BasicStructure;

  BasicStructure prim = test::ZrO_prim();

  // UnitCellCoord is a data structure representing crystal sites via integral coordinates:
  // - the sublattice index, (b, index into Structure::basis)
  // - the unit cell coordinates, (i,j,k, multiples of the unit cell vectors, stored as UnitCell object)

  UnitCellCoord uccoord_1 {0, UnitCell {1, 2, 3} };
  EXPECT_EQ(uccoord_1.sublattice(), 0);
  EXPECT_EQ(uccoord_1.unitcell()[0], 1);
  EXPECT_EQ(uccoord_1.unitcell()[1], 2);
  EXPECT_EQ(uccoord_1.unitcell()[2], 3);

  // There is also a 4-integer coordinate constructor that takes {b, i, j, k};
  UnitCellCoord uccoord_2 {0, 1, 2, 3};
  EXPECT_EQ(uccoord_1, uccoord_2);

  // UnitCellCoord can be translated by a UnitCell translation vector
  uccoord_2 += UnitCell {1, 1, 1};
  EXPECT_EQ(uccoord_2.sublattice(), 0);
  EXPECT_EQ(uccoord_2.unitcell()[0], 2);
  EXPECT_EQ(uccoord_2.unitcell()[1], 3);
  EXPECT_EQ(uccoord_2.unitcell()[2], 4);

  // UnitCellCoord can be converted to Coordinate by referring to the reference structure
  EXPECT_TRUE(uccoord_1.coordinate(prim).almost_equal(
                uccoord_1.unitcell().coordinate(prim.lattice()) + uccoord_1.sublattice_site(prim)));
}
