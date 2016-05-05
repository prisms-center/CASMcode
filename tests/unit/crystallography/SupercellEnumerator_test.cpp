#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/crystallography/Lattice.hh"

/// What is being used to test it:
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/external/Eigen/Dense"

using namespace CASM;

void autofail() {
  BOOST_CHECK_EQUAL(1, 0);
  return;
}

void hermite_init() {
  int dims = 5;
  int det = 30;

  HermiteCounter hermit_test(det, dims);

  Eigen::VectorXi init_diagonal(Eigen::VectorXi::Ones(dims));
  init_diagonal(0) = det;

  BOOST_CHECK_EQUAL(init_diagonal, hermit_test.diagonal());
  BOOST_CHECK_EQUAL(2, hermit_test.att());
  BOOST_CHECK_EQUAL(0, hermit_test.pos());


  return;
}
BOOST_AUTO_TEST_SUITE(SupercellEnumeratorTest)

BOOST_AUTO_TEST_CASE(HermitConstruction) {
  hermite_init();
}

BOOST_AUTO_TEST_SUITE_END()
