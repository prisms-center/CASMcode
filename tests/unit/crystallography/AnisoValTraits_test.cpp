#include "casm/crystallography/AnisoValTraits.hh"

#include "casm/clex/SimpleStructureTools.hh"
#include "gtest/gtest.h"

using namespace CASM;
using namespace CASM::clex_SimpleStructureTools_impl;

/// Assert to make sure that when std::set<TransfromDirective> is made with
/// atomize and an AnisoValTrait (which should be applied before atomize),
/// AnisoValTrait exists before atomize in the set
void assert_anisovaltrait_before_atomize(const std::string& anisovaltrait) {
  std::set<TransformDirective> transform_directives(
      {TransformDirective("atomize"), TransformDirective(anisovaltrait)});
  EXPECT_EQ(transform_directives.size(), 2);

  int i = 0;
  for (const auto& directive : transform_directives) {
    if (i == 0) {
      EXPECT_EQ(directive.name(), anisovaltrait);
    }

    if (i == 1) {
      EXPECT_EQ(directive.name(), "atomize");
    }
    ++i;
  }
}

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

TEST(AnisoValTraitsTest, Cmagpsin) {
  AnisoValTraits Cmagpsin = AnisoValTraits::Cmagspin();
  std::set<std::string> must_apply_after = Cmagpsin.must_apply_after();
  EXPECT_EQ(must_apply_after.size(), 1);

  for (const std::string& atomize : must_apply_after) {
    EXPECT_EQ(atomize, "atomize");
  }
  assert_anisovaltrait_before_atomize("Cmagspin");
}

TEST(AnisoValTraitsTest, Cunitmagspin) {
  AnisoValTraits Cunitmagspin = AnisoValTraits::Cunitmagspin();
  std::set<std::string> must_apply_after = Cunitmagspin.must_apply_after();
  EXPECT_EQ(must_apply_after.size(), 1);

  for (const std::string& atomize : must_apply_after) {
    EXPECT_EQ(atomize, "atomize");
  }
  assert_anisovaltrait_before_atomize("Cunitmagspin");
}

TEST(AnisoValTraitsTest, NCmagpsin) {
  AnisoValTraits NCmagpsin = AnisoValTraits::NCmagspin();
  std::set<std::string> must_apply_after = NCmagpsin.must_apply_after();
  EXPECT_EQ(must_apply_after.size(), 1);

  for (const std::string& atomize : must_apply_after) {
    EXPECT_EQ(atomize, "atomize");
  }
  assert_anisovaltrait_before_atomize("NCmagspin");
}

TEST(AnisoValTraitsTest, NCunitmagspin) {
  AnisoValTraits NCunitmagspin = AnisoValTraits::NCunitmagspin();
  std::set<std::string> must_apply_after = NCunitmagspin.must_apply_after();
  EXPECT_EQ(must_apply_after.size(), 1);

  for (const std::string& atomize : must_apply_after) {
    EXPECT_EQ(atomize, "atomize");
  }
  assert_anisovaltrait_before_atomize("NCunitmagspin");
}

TEST(AnisoValTraitsTest, SOmagspin) {
  AnisoValTraits SOmagspin = AnisoValTraits::NCunitmagspin();
  std::set<std::string> must_apply_after = SOmagspin.must_apply_after();
  EXPECT_EQ(must_apply_after.size(), 1);

  for (const std::string& atomize : must_apply_after) {
    EXPECT_EQ(atomize, "atomize");
  }
  assert_anisovaltrait_before_atomize("SOmagspin");
}

TEST(AnisoValTraitsTest, SOunitmagspin) {
  AnisoValTraits SOunitmagspin = AnisoValTraits::SOunitmagspin();
  std::set<std::string> must_apply_after = SOunitmagspin.must_apply_after();
  EXPECT_EQ(must_apply_after.size(), 1);

  for (const std::string& atomize : must_apply_after) {
    EXPECT_EQ(atomize, "atomize");
  }
  assert_anisovaltrait_before_atomize("SOunitmagspin");
}
