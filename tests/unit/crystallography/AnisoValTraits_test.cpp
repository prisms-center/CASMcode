#include "casm/crystallography/AnisoValTraits.hh"

#include "gtest/gtest.h"

using namespace CASM;

TEST(AnisoValTraitsTest, dOrbitalOccupationConstruction) {
  AnisoValTraits occ = AnisoValTraits::d_orbital_occupation();
  EXPECT_EQ(occ.name(), "dorbitaloccupation") << "Name is wrong!";
  EXPECT_EQ(occ.dim(), 15) << "Not 15-dimensional!";
}

TEST(AnisoValTraitsTest, dOrbitalOccupationSpinPolarizedConstruction) {
  AnisoValTraits occ = AnisoValTraits::d_orbital_occupation_spin_polarized();
  EXPECT_EQ(occ.name(), "dorbitaloccupationspinpolarized") << "Name is wrong!";
  EXPECT_EQ(occ.dim(), 30) << "Not 30-dimensional!";
}
