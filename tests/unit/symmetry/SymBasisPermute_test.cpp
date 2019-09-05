#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/symmetry/SymBasisPermute.hh"

/// What is being used to test it:
#include "Common.hh"
#include "TestConfiguration.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/misc/CASM_Eigen_math.hh"

using namespace CASM;

namespace {
  struct TestConfig0 : test::TestConfiguration {

    TestConfig0(const PrimClex &primclex) :
      TestConfiguration(
        primclex,
        Eigen::Vector3i(2, 1, 1).asDiagonal(), {
      0, 0,  0, 0,  1, 1,  0, 0
    }) {

      BOOST_CHECK_EQUAL(this->scel_fg().size(), 16);
      BOOST_CHECK_EQUAL(this->config_sym_fg().size(), 8);
    }

  };
}

BOOST_AUTO_TEST_SUITE(SymBasisPermuteTest)

BOOST_AUTO_TEST_CASE(Test0) {

  test::ZrOProj proj;
  proj.check_init();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);

  TestConfig0 td(primclex);

  UnitCellCoord u(primclex.prim(), 2, 0, 0, 0);

  for(Index l = 0; l < td.config.size(); ++l) {

    UnitCellCoord u = td.config.uccoord(l);

    //std::cout << "coord: " << u.coordinate().const_cart().transpose() << std::endl;
    const SymGroup &g = td.config_sym_fg();
    BOOST_CHECK_EQUAL(true, true);
    for(Index op_i = 0; op_i < g.size(); ++op_i) {

      // check that applying symmetry via UnitCellCoord and via Coordinate give same result
      BOOST_CHECK_EQUAL(
        true,
        almost_equal(
          copy_apply(g[op_i], u.coordinate()).const_cart(),
          copy_apply(g[op_i], u).coordinate().const_cart(),
          primclex.crystallography_tol()));
    }
  }

}

BOOST_AUTO_TEST_SUITE_END()

