#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/crystallography/UnitCellCoord.hh"

/// What is being used to test it:
#include "FCCTernaryProj.hh"

using namespace CASM;


BOOST_AUTO_TEST_SUITE(UnitCellCoordTest)

BOOST_AUTO_TEST_CASE(Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();
  proj.check_composition();

  PrimClex primclex(proj.dir, default_log());

  {
    UnitCellCoord uccoord(primclex.prim(), 0, -1, 1, 1);
    Eigen::Vector3d vec = uccoord.site().cart();
    BOOST_CHECK_EQUAL(almost_equal(vec, Eigen::Vector3d(4., 0., 0.)), true);
  }

}


BOOST_AUTO_TEST_SUITE_END()
