#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/ScelEnum.hh"
//#include "casm/clex/ScelEnum_impl.hh"

/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(ScelEnumTest)

BOOST_AUTO_TEST_CASE(Test1) {

  test::ZrOProj proj;
  make_project(proj);

  PrimClex primclex(proj.dir, null_log());

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  std::vector<std::string> m_names;

  // -- Test ScelEnumByProps --------------------
  {
    ScelEnumProps enum_props(1, 5);
    ScelEnumByProps e(primclex, enum_props);

    BOOST_CHECK_EQUAL(e.name(), "ScelEnumByProps");

    auto it = e.begin();
    BOOST_CHECK(true);
    BOOST_CHECK_EQUAL(it.name(), "ScelEnumByProps");

    auto end = e.end();
    BOOST_CHECK(true);

    Index count = 0;
    for(; it != end; ++it, ++count) {
      m_names.push_back(it->name());
      //std::cout << it->name() << std::endl;
    }
    BOOST_CHECK_EQUAL(count, 20);
    BOOST_CHECK(it == end);
  }

  // -- use results to Test ScelEnumByName --------------------
  {
    ScelEnumByNameT<false> e(primclex, m_names.begin(), m_names.end());
    BOOST_CHECK_EQUAL(e.name(), "ScelEnumByName");

    auto it = e.begin();
    BOOST_CHECK(true);
    BOOST_CHECK_EQUAL(it.name(), "ScelEnumByName");

    auto end = e.end();
    BOOST_CHECK(true);

    Index count = 0;
    for(; it != end; ++it, ++count) {
      //std::cout << it->name() << std::endl;
    }
    BOOST_CHECK_EQUAL(count, 20);
    BOOST_CHECK(it == end);
  }

  rm_project(proj);
}

BOOST_AUTO_TEST_SUITE_END()
