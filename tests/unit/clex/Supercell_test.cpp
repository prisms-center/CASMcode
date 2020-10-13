#include "gtest/gtest.h"

/// What is being tested:
#include "casm/clex/Supercell_impl.hh"

/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/app/AppIO.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/LatticeIsEquivalent.hh"

using namespace CASM;

TEST(SupercellTest, Constructor1) {

  // basic construction: shared prim structure and
  //                     prim lattice -> supercell lattice transformation matrix

  auto shared_prim = std::make_shared<Structure const>(test::FCC_ternary_prim());
  Eigen::Matrix3l T = Eigen::Matrix3l::Identity();
  auto shared_supercell = std::make_shared<Supercell const>(shared_prim, T);
  EXPECT_EQ(shared_supercell->transf_mat(), T);
}

TEST(SupercellTest, Constructor2) {

  // basic construction: shared prim structure and
  //                     supercell lattice

  auto shared_prim = std::make_shared<Structure const>(test::FCC_ternary_prim());
  Eigen::Matrix3l T = Eigen::Matrix3l::Identity();
  Lattice supercell_lattice = make_superlattice(shared_prim->lattice(), T);
  auto shared_supercell = std::make_shared<Supercell const>(shared_prim, supercell_lattice);
  EXPECT_EQ(shared_supercell->lattice(), supercell_lattice);
  EXPECT_EQ(shared_supercell->transf_mat(), T);
}

TEST(SupercellTest, ConstructorFail) {

  // fail construction: shared prim structure with no basis sites

  auto shared_prim = std::make_shared<Structure const>(test::no_basis_prim());
  Eigen::Matrix3l T = Eigen::Matrix3l::Identity();
  ASSERT_ANY_THROW(std::make_shared<Supercell const>(shared_prim, T));
}

TEST(SupercellTest, TestSupercellName) {

  ScopedNullLogging logging;
  test::FCCTernaryProj proj;
  proj.check_init();
  PrimClex primclex(proj.dir);

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  {
    Supercell scel {&primclex, Lattice(a, b, c)};
    EXPECT_EQ(true, true);
    EXPECT_EQ(scel.name(), "SCEL1_1_1_1_0_0_0");
  }

  {
    // standard cubic FCC unit cell
    Supercell scel {&primclex, Lattice(c + b - a, a - b + c, a + b - c)};
    EXPECT_EQ(true, true);
    EXPECT_EQ(scel.name(), "SCEL4_2_2_1_1_1_0");
  }

  {
    // non-standard, but equivalent cubic FCC unit cell
    Supercell scel {&primclex, Lattice(c + b - a, a - b + c, (a + b - c) + (c + b - a))};
    EXPECT_EQ(true, true);
    EXPECT_EQ(scel.name(), "SCEL4_2_2_1_1_1_0");
  }

  {
    // non-standard, transformed cubic FCC unit cell (standard w/ z+'c')
    Lattice lat(c + b - a, a - b + c, 2 * (a + b - c));
    Lattice rotated_lat = sym::copy_apply(primclex.prim().point_group()[10], lat);

    Supercell scel {&primclex, rotated_lat};
    EXPECT_EQ(true, true);
    EXPECT_EQ(scel.name(), "SCEL8_4_2_1_1_3_2.5");
  }

}
