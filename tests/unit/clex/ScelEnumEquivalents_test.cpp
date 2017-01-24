#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/ScelEnumEquivalents.hh"

/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(ScelEnumEquivalentsTest)

BOOST_AUTO_TEST_CASE(Test1) {

  test::ZrOProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.get_prim().lattice().vectors();

  {
    Supercell scel {&primclex, Lattice(a, b, c)};
    ScelEnumEquivalents e(scel);
    BOOST_CHECK_EQUAL(1, std::distance(e.begin(), e.end()));
  }

  {
    Supercell scel {&primclex, Lattice(2.*a, b, c)};
    ScelEnumEquivalents e(scel);
    BOOST_CHECK_EQUAL(3, std::distance(e.begin(), e.end()));
  }

  {
    Supercell scel {&primclex, Lattice(2.*a, 2.*b, c)};
    ScelEnumEquivalents e(scel);
    BOOST_CHECK_EQUAL(1, std::distance(e.begin(), e.end()));
  }

}

BOOST_AUTO_TEST_CASE(Test2) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.get_prim().lattice().vectors();

  {
    Supercell scel {&primclex, Lattice(a, b, c)};
    ScelEnumEquivalents e(scel);
    BOOST_CHECK_EQUAL(1, std::distance(e.begin(), e.end()));
  }

  {
    Supercell scel {&primclex, Lattice(2.*a, b, c)};
    ScelEnumEquivalents e(scel);
    BOOST_CHECK_EQUAL(4, std::distance(e.begin(), e.end()));
  }

  {
    Supercell scel {&primclex, Lattice(c, a - b, a + b - c)};
    ScelEnumEquivalents e(scel);
    BOOST_CHECK_EQUAL(3, std::distance(e.begin(), e.end()));
  }

  {
    Supercell scel {&primclex, Lattice(2.*a, 2.*b, c)};
    ScelEnumEquivalents e(scel);
    BOOST_CHECK_EQUAL(6, std::distance(e.begin(), e.end()));
  }

  {
    Supercell scel {&primclex, Lattice(2.*a, 2.*b, 2.*c)};
    ScelEnumEquivalents e(scel);
    BOOST_CHECK_EQUAL(1, std::distance(e.begin(), e.end()));
  }

  {
    Supercell scel {&primclex, Lattice(4.*a, 2.*b, 1.*c)};
    ScelEnumEquivalents e(scel);
    BOOST_CHECK_EQUAL(12, std::distance(e.begin(), e.end()));
  }

}

BOOST_AUTO_TEST_SUITE_END()
