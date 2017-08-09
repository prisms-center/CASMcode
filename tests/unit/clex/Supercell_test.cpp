#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/Supercell_impl.hh"

/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/app/AppIO.hh"
#include "casm/crystallography/Structure.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(SupercellTest)

BOOST_AUTO_TEST_CASE(TestSupercellName) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  {
    Supercell scel {&primclex, Lattice(a, b, c)};
    std::cout << "lat 0: \n" << scel.lattice().lat_column_mat() << std::endl;
    std::cout << "scel.name() 0: " << scel.name() << std::endl;
    BOOST_CHECK_EQUAL(true, true);
  }

  {
    // standard cubic FCC unit cell
    Supercell scel {&primclex, Lattice(c + b - a, a - b + c, a + b - c)};
    std::cout << "lat 1: \n" << scel.lattice().lat_column_mat() << std::endl;
    std::cout << "scel.name() 1: " << scel.name() << std::endl;
    BOOST_CHECK_EQUAL(true, true);
  }

  {
    // non-standard cubic FCC unit cell (standard w/ z+'c'
    Supercell scel {&primclex, Lattice(c + b - a, a - b + c, (a + b - c) + c)};
    std::cout << "lat 2: \n" << scel.lattice().lat_column_mat() << std::endl;
    std::cout << "scel.name() 2: " << scel.name() << std::endl;
    BOOST_CHECK_EQUAL(true, true);

    Index i = 0;
    for(const auto &op : scel.prim().point_group()) {
      std::cout << " op: " << i << "  index: " << op.index() << std::endl;
      ++i;
    }
    std::cout << "from_canonical: " << scel.from_canonical().index() << std::endl;

  }

}

BOOST_AUTO_TEST_SUITE_END()
