#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/Supercell_impl.hh"

/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/app/AppIO.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/LatticeIsEquivalent.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(SupercellTest)

BOOST_AUTO_TEST_CASE(TestSupercellName) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  {
    Supercell scel {&primclex, Lattice(a, b, c)};
    BOOST_CHECK_EQUAL(true, true);
    BOOST_CHECK_EQUAL(scel.name(), "SCEL1_1_1_1_0_0_0");
  }

  {
    // standard cubic FCC unit cell
    Supercell scel {&primclex, Lattice(c + b - a, a - b + c, a + b - c)};
    BOOST_CHECK_EQUAL(true, true);
    BOOST_CHECK_EQUAL(scel.name(), "SCEL4_2_2_1_1_1_0");
  }

  {
    // non-standard, but equivalent cubic FCC unit cell
    Supercell scel {&primclex, Lattice(c + b - a, a - b + c, (a + b - c) + (c + b - a))};
    BOOST_CHECK_EQUAL(true, true);
    BOOST_CHECK_EQUAL(scel.name(), "SCEL4_2_2_1_1_1_0");
  }

  {
    // non-standard, transformed cubic FCC unit cell (standard w/ z+'c')
    Lattice lat(c + b - a, a - b + c, 2 * (a + b - c));
    Lattice rotated_lat = copy_apply(primclex.prim().point_group()[10], lat);

    Supercell scel {&primclex, rotated_lat};
    BOOST_CHECK_EQUAL(true, true);
    BOOST_CHECK_EQUAL(scel.name(), "SCEL8_4_2_1_1_3_2.5");
  }

}

BOOST_AUTO_TEST_SUITE_END()
