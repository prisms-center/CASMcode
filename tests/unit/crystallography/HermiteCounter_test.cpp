#include "casm/crystallography/HermiteCounter.hh"

#include "gtest/gtest.h"

using namespace CASM;
using namespace xtal;

void hermite_init() {
  int dims = 5;
  int det = 30;

  HermiteCounter hermit_test(det, dims);

  Eigen::VectorXi init_diagonal(Eigen::VectorXi::Ones(dims));
  init_diagonal(0) = det;

  EXPECT_EQ(init_diagonal, hermit_test.diagonal());
  EXPECT_EQ(0, hermit_test.position());

  auto tricounter =
      HermiteCounter_impl::_upper_tri_counter(hermit_test.diagonal());
  Eigen::VectorXi startcount(
      Eigen::VectorXi::Zero(HermiteCounter_impl::upper_size(dims)));
  EXPECT_EQ(tricounter.current(), startcount);

  Eigen::VectorXi endcount(
      Eigen::VectorXi::Zero(HermiteCounter_impl::upper_size(dims)));
  endcount(0) = (det - 1);
  endcount(1) = (det - 1);
  endcount(2) = (det - 1);
  endcount(3) = (det - 1);

  auto finalcounterstate = tricounter;

  for (; tricounter.valid(); ++tricounter) {
    finalcounterstate = tricounter;
  }

  EXPECT_EQ(finalcounterstate.current(), endcount);

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
  EXPECT_EQ(p0, p + 1);
  EXPECT_EQ(diagonal0(p), 1);
  EXPECT_EQ(diagonal0(p + 1), 2);

  p = 3;
  diagonal1(p) = 6;
  int p1 = HermiteCounter_impl::_spill_factor(diagonal1, p, 2);
  EXPECT_EQ(p1, p + 1);
  EXPECT_EQ(diagonal1(p), 3);
  EXPECT_EQ(diagonal1(p + 1), 2);

  p = 3;
  diagonal2(p) = 6;
  int p2 = HermiteCounter_impl::_spill_factor(diagonal2, p, 4);
  EXPECT_EQ(p2, p + 1);
  EXPECT_EQ(diagonal2(p), 1);
  EXPECT_EQ(diagonal2(p + 1), 6);

  p = 2;
  diagonal3(p) = 8;
  int p3 = HermiteCounter_impl::_spill_factor(diagonal3, p, 4);
  EXPECT_EQ(p3, p + 1);
  EXPECT_EQ(diagonal3(p), 2);
  EXPECT_EQ(diagonal3(p + 1), 4);

  return;
}

void next_position_test() {
  // Example increment from one possible diagonal to the next
  Eigen::VectorXi diagonal(Eigen::VectorXi::Ones(5));
  Eigen::VectorXi next_diagonal(Eigen::VectorXi::Ones(5));
  diagonal(0) = 6;
  next_diagonal(0) = 3;
  next_diagonal(1) = 2;
  int p = 0;

  p = HermiteCounter_impl::next_spill_position(diagonal, p);

  EXPECT_EQ(diagonal, next_diagonal);
  EXPECT_EQ(p, 1);

  diagonal = Eigen::VectorXi::Ones(5);
  next_diagonal = Eigen::VectorXi::Ones(5);
  //[1 2 1 1 3]
  diagonal(1) = 2;
  diagonal(4) = 3;
  //[1 1 6 1 1]
  next_diagonal(2) = 6;

  p = 4;
  p = HermiteCounter_impl::next_spill_position(diagonal, p);

  EXPECT_EQ(diagonal, next_diagonal);
  EXPECT_EQ(p, 2);

  //*************/
  // Make sure every enumerated diagonal has the right determinant

  int det = 2 * 3 * 5 * 7;
  int dims = 5;

  Eigen::VectorXi diag = Eigen::VectorXi::Ones(dims);
  diag(0) = det;

  p = 0;
  while (p != diag.size()) {
    int testdet = 1;
    for (int i = 0; i < diag.size(); i++) {
      testdet = testdet * diag(i);
    }
    EXPECT_EQ(det, testdet);
    p = CASM::xtal::HermiteCounter_impl::next_spill_position(diag, p);
  }

  return;
}

void triangle_count_test() {
  HermiteCounter::Index totals = HermiteCounter_impl::upper_size(7);
  EXPECT_EQ(totals, -7 + 7 + 6 + 5 + 4 + 3 + 2 + 1);

  int dims = 5;
  // int det = 30;

  Eigen::VectorXi mid_diagonal(Eigen::VectorXi::Ones(dims));
  mid_diagonal(0) = 5;
  mid_diagonal(1) = 3;
  mid_diagonal(4) = 2;

  auto countertest = HermiteCounter_impl::_upper_tri_counter(mid_diagonal);
  auto finalcount = countertest;

  for (; countertest.valid(); countertest++) {
    finalcount = countertest;
  }

  // The initial matrix is 5x5 with diagonal [ 5 3 1 1 2 ], so it has
  // determinant=30 The Hermite matrix with highest ranking for this determinant
  // will therefore be:
  //    5 4 4 4 4
  //    0 3 2 2 2
  //    0 0 1 0 0
  //    0 0 0 1 0
  //    0 0 0 0 2
  // Which gives the upper triangular vector [ 4 4 4 4 2 2 2 0 0 0 ]

  Eigen::VectorXi end_count_value(
      Eigen::VectorXi::Zero(HermiteCounter_impl::upper_size(5)));
  end_count_value(0) = 4;
  end_count_value(1) = 4;
  end_count_value(2) = 4;
  end_count_value(3) = 4;
  end_count_value(4) = 2;
  end_count_value(5) = 2;
  end_count_value(6) = 2;
  // Rest of the values are zero

  EXPECT_EQ(finalcount.current(), end_count_value);

  return;
}

void matrix_construction_test() {
  Eigen::VectorXi diag;
  diag.resize(4);
  diag << 2, 4, 6, 8;
  Eigen::VectorXi upper;
  upper.resize(3 + 2 + 1);
  upper << 11, 12, 13, 21, 22, 33;

  Eigen::MatrixXi diagmat;
  diagmat.resize(4, 4);
  diagmat << 2, 11, 12, 13, 0, 4, 21, 22, 0, 0, 6, 33, 0, 0, 0, 8;

  EXPECT_EQ(diagmat, HermiteCounter_impl::_zip_matrix(diag, upper));

  return;
}

void increment_test() {
  HermiteCounter hermit_test(6, 4);

  Eigen::MatrixXi hermmat;
  hermmat.resize(4, 4);

  // Test starting status
  hermmat << 6, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;

  EXPECT_EQ(hermmat, hermit_test());

  // Test next status
  ++hermit_test;
  hermmat << 6, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;

  EXPECT_EQ(hermmat, hermit_test());

  // Jump to just before you need a new diagonal
  hermmat << 6, 5, 5, 5, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;

  while (hermit_test() != hermmat) {
    ++hermit_test;
  }

  // Check diagonal jump
  ++hermit_test;
  hermmat << 3, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;

  EXPECT_EQ(hermmat, hermit_test());

  // Check invalidation and last status
  auto lastherm = hermmat;
  while (hermit_test.determinant() != 7) {
    lastherm = hermit_test();
    ++hermit_test;
  }

  hermmat << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 6;

  EXPECT_EQ(hermmat, lastherm);

  // Check determinant jump
  hermit_test = HermiteCounter(3, 4);

  // Jump to just before you need a new determinant

  hermmat << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 3;

  while (hermit_test() != hermmat) {
    ++hermit_test;
  }

  // Check determinant jump
  ++hermit_test;

  hermmat << 4, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;

  EXPECT_EQ(hermmat, hermit_test());
  return;
}

void reset_test() {
  HermiteCounter hermit_test(1, 3);

  Eigen::MatrixXi hermmat;
  hermmat.resize(3, 3);

  Eigen::MatrixXi startmat = hermit_test();

  // Skip to one of the bigger determinants
  hermmat << 2, 1, 1, 0, 2, 1, 0, 0, 1;

  while (hermit_test() != hermmat) {
    ++hermit_test;
  }

  hermmat << 4, 0, 0, 0, 1, 0, 0, 0, 1;

  hermit_test.reset_current();

  EXPECT_EQ(hermmat, hermit_test());

  hermit_test.jump_to_determinant(1);

  EXPECT_EQ(startmat, hermit_test());

  return;
}

void expand_dims_test() {
  Eigen::MatrixXi expandmat(Eigen::MatrixXi::Ones(5, 5));
  expandmat = expandmat * 3;

  Eigen::VectorXi expanddims(8);
  expanddims << 1, 1, 1, 0, 1, 0, 0, 1;

  Eigen::MatrixXi expandedmat(8, 8);
  expandmat << 3, 3, 3, 0, 3, 0, 0, 3, 3, 3, 3, 0, 3, 0, 0, 3, 3, 3, 3, 0, 3, 0,
      0, 3, 0, 0, 0, 1, 0, 0, 0, 0, 3, 3, 3, 0, 3, 0, 0, 3, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0, 3, 3, 3, 0, 3, 0, 0, 1;

  EXPECT_EQ(expandmat,
            HermiteCounter_impl::_expand_dims_old(expandmat, expanddims));

  HermiteCounter minicount(1, 4);
  for (int i = 0; i < 12; i++) {
    ++minicount;
  }

  Eigen::MatrixXi endcount(4, 4);
  endcount << 1, 0, 0, 0, 0, 2, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1;

  EXPECT_EQ(endcount, minicount());

  Eigen::MatrixXi transmat(Eigen::MatrixXi::Identity(6, 6));

  Eigen::MatrixXi expanded =
      HermiteCounter_impl::_expand_dims(minicount(), transmat);
  Eigen::MatrixXi blockmat(6, 6);
  blockmat << 1, 0, 0, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1;

  EXPECT_EQ(blockmat, expanded);

  Eigen::Matrix2Xi miniherm;
  miniherm << 2, 1, 0, 3;

  Eigen::Matrix3i minitrans;
  minitrans << 1, 0, 0, 0, 0, 1, 0, 1, 0;

  Eigen::Matrix3i miniexpand;
  miniexpand << 2, 1, 0, 0, 0, 1, 0, 3, 0;

  EXPECT_EQ(HermiteCounter_impl::_expand_dims(miniherm, minitrans), miniexpand);

  return;
}

void unroll_test() {
  Eigen::MatrixXi mat5(5, 5);
  mat5 << 1, 12, 11, 10, 9, 0, 2, 13, 15, 8, 0, 0, 3, 14, 7, 0, 0, 0, 4, 6, 0,
      0, 0, 0, 5;

  auto result5 = HermiteCounter_impl::_canonical_unroll(mat5);

  for (int i = 0; i < result5.size(); ++i) {
    EXPECT_EQ(i + 1, result5(i));
  }

  /* Eigen::MatrixXi mat3(3, 3); */
  Eigen::Matrix3i mat3;
  mat3 << 1, 6, 5, 0, 2, 4, 0, 0, 3;

  auto result3 = HermiteCounter_impl::_canonical_unroll(mat3);
  for (int i = 0; i < result3.size(); ++i) {
    std::cerr << result3(i) << std::endl;
    EXPECT_EQ(i + 1, result5(i));
  }

  return;
}

void compare_test() {
  Eigen::Matrix3i low, high;

  low << 1, 9, 9, 0, 9, 9, 0, 9, 9;

  high << 2, 0, 0, 0, 1, 0, 0, 0, 1;

  EXPECT_TRUE(HermiteCounter_impl::_canonical_compare(low, high));

  low << 1, 9, 9, 0, 9, 9, 0, 9, 9;

  high << 1, 10, 9, 0, 9, 9, 0, 9, 9;

  EXPECT_TRUE(HermiteCounter_impl::_canonical_compare(low, high));

  return;
}

TEST(HermiteCounterTest, HermiteConstruction) { hermite_init(); }

TEST(HermiteCounterTest, HermiteImpl_spill) { spill_test(); }

TEST(HermiteCounterTest, HermiteImpl_next_position) { next_position_test(); }

TEST(HermiteCounterTest, HermiteImpl_triangle_count) { triangle_count_test(); }

TEST(HermiteCounterTest, HermiteImpl_matrix_construction) {
  matrix_construction_test();
}

TEST(HermiteCounterTest, HermiteImpl_reset_test) { reset_test(); }

TEST(HermiteCounterTest, HermiteImpl_unroll_test) { unroll_test(); }

TEST(HermiteCounterTest, HermiteImpl_compare_test) { compare_test(); }

TEST(HermiteCounterTest, HermiteCounting_increment_test) { increment_test(); }
