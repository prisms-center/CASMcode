#include "gtest/gtest.h"

/// What is being tested:
#include "casm/crystallography/LatticeEnumEquivalents.hh"

/// What is being used to test it:
#include "casm/crystallography/Lattice_impl.hh"
#include "casm/crystallography/Structure.hh"
#include "ZrOProj.hh"

using namespace CASM;
using namespace test;

TEST(LatticeEnumEquivalentsTest, Test1) {

  Structure ZrO(ZrO_prim());

  LatticeEnumEquivalents enumerator(ZrO.lattice(), ZrO.factor_group());
  EXPECT_TRUE(1) << "LatticeEnumEquivalents construction failed";

  auto begin = enumerator.begin();
  EXPECT_TRUE(1) << "LatticeEnumEquivalents::begin() failed";

  auto end = enumerator.end();
  EXPECT_TRUE(1) << "LatticeEnumEquivalents::end() failed";

  // for prim, should only be one equivalent
  EXPECT_EQ(1, std::distance(begin, end));

  // LatticeEnumEquivalents is an InputEnum and allows only a single pass
  EXPECT_EQ(enumerator.valid(), false);

  // LatticeEnumEquivalents is an InputEnum and allows only a single pass
  EXPECT_EQ(0, std::distance(enumerator.begin(), enumerator.end()));

}

TEST(LatticeEnumEquivalentsTest, Test2) {

  Structure ZrO(ZrO_prim());

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = ZrO.lattice().vectors();

  {
    LatticeEnumEquivalents e(Lattice(2.*a, b, c), ZrO.factor_group());
    EXPECT_EQ(3, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(2.*a, 2 * b, c), ZrO.factor_group());
    EXPECT_EQ(1, std::distance(e.begin(), e.end()));
  }

}

TEST(LatticeEnumEquivalentsTest, Test3) {

  Lattice lat = Lattice::hexagonal();

  SymGroup pg = SymGroup::lattice_point_group(lat);
  //lat.generate_point_group(pg);

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = lat.vectors();

  {
    LatticeEnumEquivalents e(lat, pg);
    EXPECT_EQ(1, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(2.*a, b, c), pg);
    EXPECT_EQ(3, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(2.*a, 2 * b, c), pg);
    EXPECT_EQ(1, std::distance(e.begin(), e.end()));
  }

}

TEST(LatticeEnumEquivalentsTest, Test4) {

  Lattice lat = Lattice::cubic();
  SymGroup pg = SymGroup::lattice_point_group(lat);

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = lat.vectors();

  {
    LatticeEnumEquivalents e(lat, pg);
    EXPECT_EQ(1, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(2.*a, b, c), pg);
    EXPECT_EQ(3, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(2.*a, 2.*b, c), pg);
    EXPECT_EQ(3, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(2.*a, 2.*b, 2.*c), pg);
    EXPECT_EQ(1, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(4.*a, 2.*b, 1.*c), pg);
    EXPECT_EQ(6, std::distance(e.begin(), e.end()));
  }

}

TEST(LatticeEnumEquivalentsTest, Test5) {

  Lattice lat = Lattice::fcc();

  SymGroup pg = SymGroup::lattice_point_group(lat);

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = lat.vectors();

  {
    LatticeEnumEquivalents e(lat, pg);
    EXPECT_EQ(1, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(2.*a, b, c), pg);
    EXPECT_EQ(4, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(c, a - b, a + b - c), pg);
    EXPECT_EQ(3, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(2.*a, 2.*b, c), pg);
    EXPECT_EQ(6, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(2.*a, 2.*b, 2.*c), pg);
    EXPECT_EQ(1, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(4.*a, 2.*b, 1.*c), pg);
    EXPECT_EQ(12, std::distance(e.begin(), e.end()));
  }

}

