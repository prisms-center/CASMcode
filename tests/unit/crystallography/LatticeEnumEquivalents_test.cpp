#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/crystallography/LatticeEnumEquivalents.hh"

/// What is being used to test it:
#include "casm/crystallography/Lattice_impl.hh"
#include "casm/crystallography/Structure.hh"
#include "ZrOProj.hh"

using namespace CASM;
using namespace test;

BOOST_AUTO_TEST_SUITE(LatticeEnumEquivalentsTest)

BOOST_AUTO_TEST_CASE(Test1) {

  Structure ZrO(ZrO_prim());

  LatticeEnumEquivalents enumerator(ZrO.lattice(), ZrO.factor_group());
  BOOST_CHECK_MESSAGE(1, "LatticeEnumEquivalents construction failed");

  auto begin = enumerator.begin();
  BOOST_CHECK_MESSAGE(1, "LatticeEnumEquivalents::begin() failed");

  auto end = enumerator.end();
  BOOST_CHECK_MESSAGE(1, "LatticeEnumEquivalents::end() failed");

  // for prim, should only be one equivalent
  BOOST_CHECK_EQUAL(1, std::distance(begin, end));

  // LatticeEnumEquivalents is an InputEnum and allows only a single pass
  BOOST_CHECK_EQUAL(enumerator.valid(), false);

  // LatticeEnumEquivalents is an InputEnum and allows only a single pass
  BOOST_CHECK_EQUAL(0, std::distance(enumerator.begin(), enumerator.end()));

}

BOOST_AUTO_TEST_CASE(Test2) {

  Structure ZrO(ZrO_prim());

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = ZrO.lattice().vectors();

  {
    LatticeEnumEquivalents e(Lattice(2.*a, b, c), ZrO.factor_group());
    BOOST_CHECK_EQUAL(3, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(2.*a, 2 * b, c), ZrO.factor_group());
    BOOST_CHECK_EQUAL(1, std::distance(e.begin(), e.end()));
  }

}

BOOST_AUTO_TEST_CASE(Test3) {

  Lattice lat = Lattice::hexagonal();

  MasterSymGroup pg;
  lat.generate_point_group(pg);

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = lat.vectors();

  {
    LatticeEnumEquivalents e(lat, pg);
    BOOST_CHECK_EQUAL(1, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(2.*a, b, c), pg);
    BOOST_CHECK_EQUAL(3, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(2.*a, 2 * b, c), pg);
    BOOST_CHECK_EQUAL(1, std::distance(e.begin(), e.end()));
  }

}

BOOST_AUTO_TEST_CASE(Test4) {

  Lattice lat = Lattice::cubic();

  MasterSymGroup pg;
  lat.generate_point_group(pg);

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = lat.vectors();

  {
    LatticeEnumEquivalents e(lat, pg);
    BOOST_CHECK_EQUAL(1, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(2.*a, b, c), pg);
    BOOST_CHECK_EQUAL(3, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(2.*a, 2.*b, c), pg);
    BOOST_CHECK_EQUAL(3, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(2.*a, 2.*b, 2.*c), pg);
    BOOST_CHECK_EQUAL(1, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(4.*a, 2.*b, 1.*c), pg);
    BOOST_CHECK_EQUAL(6, std::distance(e.begin(), e.end()));
  }

}

BOOST_AUTO_TEST_CASE(Test5) {

  Lattice lat = Lattice::fcc();

  MasterSymGroup pg;
  lat.generate_point_group(pg);

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = lat.vectors();

  {
    LatticeEnumEquivalents e(lat, pg);
    BOOST_CHECK_EQUAL(1, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(2.*a, b, c), pg);
    BOOST_CHECK_EQUAL(4, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(c, a - b, a + b - c), pg);
    BOOST_CHECK_EQUAL(3, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(2.*a, 2.*b, c), pg);
    BOOST_CHECK_EQUAL(6, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(2.*a, 2.*b, 2.*c), pg);
    BOOST_CHECK_EQUAL(1, std::distance(e.begin(), e.end()));
  }

  {
    LatticeEnumEquivalents e(Lattice(4.*a, 2.*b, 1.*c), pg);
    BOOST_CHECK_EQUAL(12, std::distance(e.begin(), e.end()));
  }

}

BOOST_AUTO_TEST_SUITE_END()
