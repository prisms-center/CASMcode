#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/crystallography/Site.hh"

/// What is being used to test it:
#include "casm/misc/CASM_Eigen_math.hh"

using namespace CASM;


BOOST_AUTO_TEST_SUITE(SiteTest)

BOOST_AUTO_TEST_CASE(Test1) {

  Eigen::Vector3d vec(0.0, 0.2, 0.4);
  double tol(1e-5);
  bool divisible(true);

  Eigen::Matrix3d L;
  L << 1.0, 0.0, 0.0,
  0.0, 1.0, 0.0,
  0.0, 0.0, 1.0;

  Lattice lat(L);
  Coordinate coord(vec, lat, CART);

  Site site_a(lat);
  BOOST_CHECK_EQUAL(site_a.site_occupant().size(), 0);

  Site site_b(coord, "A");
  BOOST_CHECK_EQUAL(site_b.site_occupant().size(), 1);

  Site site_c(coord, {Molecule::make_atom("A"), Molecule::make_atom("B"), Molecule::make_atom("C")});
  BOOST_CHECK_EQUAL(site_c.site_occupant().size(), 3);

  std::vector<Molecule> tocc {Molecule::make_atom("C"), Molecule::make_atom("D")};
  site_c.set_allowed_species(tocc);
  BOOST_CHECK_EQUAL(site_c.site_occupant().size(), 2);

}


BOOST_AUTO_TEST_SUITE_END()
