#include "gtest/gtest.h"
#include "casm/crystallography/AnisoValTraits.hh"

using namespace CASM;

TEST(AnisoValTraitsTest, dOrbitalOccupationConstruction) {
    AnisoValTraits occ = AnisoValTraits::d_orbital_occupation();
    EXPECT_EQ(occ.name(), "d_orbital_occupation") << "Name is wrong!";
    EXPECT_EQ(occ.dim(), 30) << "Not 30-dimensional!";
}
