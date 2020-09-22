#include "gtest/gtest.h"

/// What is being tested:
#include "casm/clex/ScelEnum.hh"

/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "ZrOProj.hh"
#include "casm/crystallography/Structure.hh"

using namespace CASM;

TEST(EnumeratorTest, Test1) {

  test::ZrOProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  std::vector<std::string> m_names;

  // -- test an input enumerator --------------------
  {
    xtal::ScelEnumProps enum_props(1, 5);
    ScelEnumByProps e(primclex, enum_props);

    EXPECT_EQ(e.name(), "ScelEnumByProps");

    auto it = e.begin();
    auto end = e.end();

    Index count = 0;
    for(; it != end; ++it, ++count) {
      m_names.push_back(it->name());
      //std::cout << it->name() << std::endl;
    }
    EXPECT_EQ(count, 20);
    EXPECT_TRUE(it == end);
    EXPECT_TRUE(!e.valid());
  }

  // -- test a random access enumerator --------------------
  ScelEnumByName e(primclex, m_names.begin(), m_names.end());
  {
    auto it = e.begin();
    EXPECT_EQ(it.step(), 0);
    ++it;
    EXPECT_EQ(it.step(), 1);

    auto it_B = e.begin();
    EXPECT_EQ(it_B.step(), 0);
    ++it_B;
    EXPECT_EQ(it_B.step(), 1);

    EXPECT_EQ(it.step(), 1);

    it += 5;
    EXPECT_EQ(it.step(), 6);

    it -= 1;
    EXPECT_EQ(it.step(), 5);
  }

  {
    EXPECT_TRUE(e.end() == e.end());
  }

  {
    auto it = e.begin();
    it += 20;
    EXPECT_TRUE(it == e.end());
  }

  {
    auto it = e.begin();
    auto it_B = e.begin();
    EXPECT_TRUE(it == it_B);

    ++it;
    EXPECT_TRUE(it != it_B);

    ++it_B;
    EXPECT_TRUE(it == it_B);

    EXPECT_EQ(std::distance(it, it_B), 0);

    it_B += 2;
    EXPECT_EQ(std::distance(it, it_B), 2);
    EXPECT_EQ(std::distance(it_B, it), -2);
  }

  {
    EXPECT_EQ(std::distance(e.begin(), e.begin()), 0);
    EXPECT_EQ(std::distance(e.end(), e.end()), 0);
    EXPECT_EQ(std::distance(e.begin(), e.end()), e.size());
    EXPECT_EQ(std::distance(e.begin(), e.end()), e.size());
  }

  {
    EXPECT_EQ(std::distance(e.rbegin(), e.rbegin()), 0);
    EXPECT_EQ(std::distance(e.rend(), e.rend()), 0);
    EXPECT_EQ(std::distance(e.rbegin(), e.rend()), e.size());
    EXPECT_EQ(std::distance(e.rbegin(), e.rend()), e.size());
  }

}
