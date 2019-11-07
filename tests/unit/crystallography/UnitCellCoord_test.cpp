#include "gtest/gtest.h"

/// What is being tested:
#include "casm/crystallography/UnitCellCoord.hh"

/// What is being used to test it:
#include "casm/misc/CASM_Eigen_math.hh"
#include "FCCTernaryProj.hh"

using namespace CASM;


TEST(UnitCellCoordTest, Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();
  proj.check_composition();

  PrimClex primclex(proj.dir, default_log());

  {
    UnitCellCoord uccoord(primclex.prim(), 0, -1, 1, 1);
    Eigen::Vector3d vec = uccoord.site().cart();
    EXPECT_EQ(almost_equal(vec, Eigen::Vector3d(4., 0., 0.)), true);
  }

}
