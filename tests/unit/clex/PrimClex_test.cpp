#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/PrimClex.hh"

/// What is being used to test it:

#include "casm/app/ProjectBuilder.hh"
#include "casm/crystallography/Structure.hh"
#include "Common.hh"
#include "FCCTernaryProj.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(PrimClexTest)

BOOST_AUTO_TEST_CASE(Basics) {

  Structure prim(test::FCC_ternary_prim());
  BOOST_CHECK_EQUAL(prim.basis().size(), 1);

  // Construct from prim
  PrimClex primclex(prim, null_log());
  BOOST_CHECK_EQUAL(primclex.prim().basis().size(), 1);

}

BOOST_AUTO_TEST_SUITE_END()
