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
  BOOST_CHECK_EQUAL(0, hermit_test.pos());

  auto tricounter = HermiteCounter_impl::_upper_tri_counter(hermit_test.diagonal());
  Eigen::VectorXi startcount(Eigen::VectorXi::Zero(HermiteCounter_impl::upper_size(dims)));
  BOOST_CHECK_EQUAL(tricounter.current(), startcount);


  Eigen::VectorXi endcount(Eigen::VectorXi::Zero(HermiteCounter_impl::upper_size(dims)));
  endcount(0) = (det - 1);
  endcount(1) = (det - 1);
  endcount(2) = (det - 1);
  endcount(3) = (det - 1);

  auto finalcounterstate = tricounter;

  for(; tricounter.valid(); ++tricounter) {
    finalcounterstate = tricounter;
  }

  BOOST_CHECK_EQUAL(finalcounterstate.current(), endcount);


  return;
}

void spill_test() {
  Eigen::VectorXi diagonal0(Eigen::VectorXi::Ones(5));
  Eigen::VectorXi diagonal1(Eigen::VectorXi::Ones(5));
  Eigen::VectorXi diagonal2(Eigen::VectorXi::Ones(5));
  Eigen::VectorXi diagonal3(Eigen::VectorXi::Ones(5));


  int p = 0;
  diagonal0(p) = 2;
  int p0 = HermiteCounter_impl::_spill_factor(diagonal0, p, 2);
  BOOST_CHECK_EQUAL(p0, p + 1);
  BOOST_CHECK_EQUAL(diagonal0(p), 1);
  BOOST_CHECK_EQUAL(diagonal0(p + 1), 2);

  p = 3;
  diagonal1(p) = 6;
  int p1 = HermiteCounter_impl::_spill_factor(diagonal1, p, 2);
  BOOST_CHECK_EQUAL(p1, p + 1);
  BOOST_CHECK_EQUAL(diagonal1(p), 3);
  BOOST_CHECK_EQUAL(diagonal1(p + 1), 2);

  p = 3;
  diagonal2(p) = 6;
  int p2 = HermiteCounter_impl::_spill_factor(diagonal2, p, 4);
  BOOST_CHECK_EQUAL(p2, p + 1);
  BOOST_CHECK_EQUAL(diagonal2(p), 1);
  BOOST_CHECK_EQUAL(diagonal2(p + 1), 6);

  p = 2;
  diagonal3(p) = 8;
  int p3 = HermiteCounter_impl::_spill_factor(diagonal3, p, 4);
  BOOST_CHECK_EQUAL(p3, p + 1);
  BOOST_CHECK_EQUAL(diagonal3(p), 2);
  BOOST_CHECK_EQUAL(diagonal3(p + 1), 4);

  return;
}

void next_position_test() {
  //Example increment from one possible diagonal to the next
  Eigen::VectorXi diagonal(Eigen::VectorXi::Ones(5));
  Eigen::VectorXi next_diagonal(Eigen::VectorXi::Ones(5));
  diagonal(0) = 6;
  next_diagonal(0) = 3;
  next_diagonal(1) = 2;
  int p = 0;

  p = HermiteCounter_impl::next_spill_position(diagonal, p);

  BOOST_CHECK_EQUAL(diagonal, next_diagonal);
  BOOST_CHECK_EQUAL(p, 1);


  diagonal = Eigen::VectorXi::Ones(5);
  next_diagonal = Eigen::VectorXi::Ones(5);
  //[1 2 1 1 3]
  diagonal(1) = 2;
  diagonal(4) = 3;
  //[1 1 6 1 1]
  next_diagonal(2) = 6;

  p = 4;
  p = HermiteCounter_impl::next_spill_position(diagonal, p);

  BOOST_CHECK_EQUAL(diagonal, next_diagonal);
  BOOST_CHECK_EQUAL(p, 2);

  //*************/
  //Make sure every enumerated diagonal has the right determinant

  int det = 2 * 3 * 5 * 7;
  int dims = 5;

  Eigen::VectorXi diag = Eigen::VectorXi::Ones(dims);
  diag(0) = det;

  p = 0;
  while(p != diag.size()) {
    int testdet = 1;
    for(int i = 0; i < diag.size(); i++) {
      testdet = testdet * diag(i);
    }
    BOOST_CHECK_EQUAL(det, testdet);
    p = CASM::HermiteCounter_impl::next_spill_position(diag, p);
  }


  return;
}

void triangle_count_test() {
  HermiteCounter::Index totals = HermiteCounter_impl::upper_size(7);
  BOOST_CHECK_EQUAL(totals, -7 + 7 + 6 + 5 + 4 + 3 + 2 + 1);

  int dims = 5;
  int det = 30;

  Eigen::VectorXi mid_diagonal(Eigen::VectorXi::Ones(dims));
  mid_diagonal(0) = 5;
  mid_diagonal(1) = 3;
  mid_diagonal(4) = 2;

  auto countertest = HermiteCounter_impl::_upper_tri_counter(mid_diagonal);
  auto finalcount = countertest;

  for(; countertest.valid(); countertest++) {
    finalcount = countertest;
  }

  Eigen::VectorXi end_count_value(Eigen::VectorXi::Zero(dims));
  end_count_value(0) = 4;
  end_count_value(1) = 4;
  end_count_value(2) = 4;
  end_count_value(3) = 4;
  end_count_value(4) = 2;
  end_count_value(5) = 2;
  end_count_value(6) = 2;

  return;
}

BOOST_AUTO_TEST_SUITE(SupercellEnumeratorTest)

BOOST_AUTO_TEST_CASE(HermiteConstruction) {
  hermite_init();
}

BOOST_AUTO_TEST_CASE(HermiteImpl) {
  spill_test();
  next_position_test();
  triangle_count_test();
}

BOOST_AUTO_TEST_SUITE_END()
