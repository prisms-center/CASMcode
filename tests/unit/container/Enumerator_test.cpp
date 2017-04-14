#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/ScelEnum.hh"
//#include "casm/clex/ScelEnum_impl.hh"

/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(EnumeratorTest)

BOOST_AUTO_TEST_CASE(Test1) {

  test::ZrOProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.get_prim().lattice().vectors();

  std::vector<std::string> m_names;

  // -- test an input enumerator --------------------
  {
    ScelEnumProps enum_props(1, 5);
    ScelEnumByProps e(primclex, enum_props);

    BOOST_CHECK_EQUAL(e.name(), "ScelEnumByProps");

    auto it = e.begin();
    auto end = e.end();

    Index count = 0;
    for(; it != end; ++it, ++count) {
      m_names.push_back(it->get_name());
      //std::cout << it->get_name() << std::endl;
    }
    BOOST_CHECK_EQUAL(count, 20);
    BOOST_CHECK(it == end);
    BOOST_CHECK(!e.valid());
  }

  // -- test a random access enumerator --------------------
  ScelEnumByNameT<false> e(primclex, m_names.begin(), m_names.end());
  {
    auto it = e.begin();
    BOOST_CHECK_EQUAL(it.step(), 0);
    ++it;
    BOOST_CHECK_EQUAL(it.step(), 1);

    auto it_B = e.begin();
    BOOST_CHECK_EQUAL(it_B.step(), 0);
    ++it_B;
    BOOST_CHECK_EQUAL(it_B.step(), 1);

    BOOST_CHECK_EQUAL(it.step(), 1);

    it += 5;
    BOOST_CHECK_EQUAL(it.step(), 6);

    it -= 1;
    BOOST_CHECK_EQUAL(it.step(), 5);
  }

  {
    BOOST_CHECK(e.end() == e.end());
  }

  {
    auto it = e.begin();
    it += 20;
    BOOST_CHECK(it == e.end());
  }

  {
    auto it = e.begin();
    auto it_B = e.begin();
    BOOST_CHECK(it == it_B);

    ++it;
    BOOST_CHECK(it != it_B);

    ++it_B;
    BOOST_CHECK(it == it_B);

    BOOST_CHECK_EQUAL(std::distance(it, it_B), 0);

    it_B += 2;
    BOOST_CHECK_EQUAL(std::distance(it, it_B), 2);
    BOOST_CHECK_EQUAL(std::distance(it_B, it), -2);
  }

  {
    BOOST_CHECK_EQUAL(std::distance(e.begin(), e.begin()), 0);
    BOOST_CHECK_EQUAL(std::distance(e.end(), e.end()), 0);
    BOOST_CHECK_EQUAL(std::distance(e.begin(), e.end()), e.size());
    BOOST_CHECK_EQUAL(std::distance(e.begin(), e.end()), e.size());
  }

  {
    BOOST_CHECK_EQUAL(std::distance(e.rbegin(), e.rbegin()), 0);
    BOOST_CHECK_EQUAL(std::distance(e.rend(), e.rend()), 0);
    BOOST_CHECK_EQUAL(std::distance(e.rbegin(), e.rend()), e.size());
    BOOST_CHECK_EQUAL(std::distance(e.rbegin(), e.rend()), e.size());
  }

}

BOOST_AUTO_TEST_SUITE_END()
