#include "gtest/gtest.h"

/// What is being tested:
#include "casm/crystallography/Site.hh"

/// What is being used to test it:
#include "casm/misc/CASM_Eigen_math.hh"

using namespace CASM;
using xtal::Coordinate;
using xtal::Lattice;
using xtal::Molecule;
using xtal::Site;

TEST(SiteTest, Test1) {
  Eigen::Vector3d vec(0.0, 0.2, 0.4);
  double tol(1e-5);
  bool divisible(true);

  Eigen::Matrix3d L;
  L << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;

  Lattice lat(L);
  Coordinate coord(vec, lat, CART);

  Site site_a(lat);
  EXPECT_EQ(site_a.occupant_dof().size(), 0);

  Site site_b(coord, "A");
  EXPECT_EQ(site_b.occupant_dof().size(), 1);

  Site site_c(coord, {Molecule::make_atom("A"), Molecule::make_atom("B"),
                      Molecule::make_atom("C")});
  EXPECT_EQ(site_c.occupant_dof().size(), 3);

  std::vector<Molecule> tocc{Molecule::make_atom("C"),
                             Molecule::make_atom("D")};
  site_c.set_allowed_occupants(tocc);
  EXPECT_EQ(site_c.occupant_dof().size(), 2);
}
