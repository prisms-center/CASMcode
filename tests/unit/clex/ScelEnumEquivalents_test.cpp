#include "gtest/gtest.h"

/// What is being tested:
#include "casm/clex/ScelEnumEquivalents.hh"

/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "ZrOProj.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/crystallography/Structure.hh"

using namespace CASM;

TEST(ScelEnumEquivalentsTest, Test1) {
  test::ZrOProj proj;
  proj.check_init();

  ScopedNullLogging logging;
  PrimClex primclex(proj.dir);

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  {
    Supercell scel{&primclex, Lattice(a, b, c)};
    ScelEnumEquivalents e(scel);
    EXPECT_EQ(1, std::distance(e.begin(), e.end()));
  }

  {
    Supercell scel{&primclex, Lattice(2. * a, b, c)};
    ScelEnumEquivalents e(scel);
    EXPECT_EQ(3, std::distance(e.begin(), e.end()));
  }

  {
    Supercell scel{&primclex, Lattice(2. * a, 2. * b, c)};
    ScelEnumEquivalents e(scel);
    EXPECT_EQ(1, std::distance(e.begin(), e.end()));
  }
}

TEST(ScelEnumEquivalentsTest, Test2) {
  test::FCCTernaryProj proj;
  proj.check_init();

  ScopedNullLogging logging;
  PrimClex primclex(proj.dir);

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  {
    Supercell scel{&primclex, Lattice(a, b, c)};
    ScelEnumEquivalents e(scel);
    EXPECT_EQ(1, std::distance(e.begin(), e.end()));
  }

  {
    Supercell scel{&primclex, Lattice(2. * a, b, c)};
    ScelEnumEquivalents e(scel);
    EXPECT_EQ(4, std::distance(e.begin(), e.end()));
  }

  {
    Supercell scel{&primclex, Lattice(c, a - b, a + b - c)};
    ScelEnumEquivalents e(scel);
    EXPECT_EQ(3, std::distance(e.begin(), e.end()));
  }

  {
    Supercell scel{&primclex, Lattice(2. * a, 2. * b, c)};
    ScelEnumEquivalents e(scel);
    EXPECT_EQ(6, std::distance(e.begin(), e.end()));
  }

  {
    Supercell scel{&primclex, Lattice(2. * a, 2. * b, 2. * c)};
    ScelEnumEquivalents e(scel);
    EXPECT_EQ(1, std::distance(e.begin(), e.end()));
  }

  {
    Supercell scel{&primclex, Lattice(4. * a, 2. * b, 1. * c)};
    ScelEnumEquivalents e(scel);
    EXPECT_EQ(12, std::distance(e.begin(), e.end()));
  }
}
