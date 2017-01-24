#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/crystallography/Lattice.hh"

/// What is being used to test it:
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/symmetry/SymGroup.hh"

using namespace CASM;

void lattice_pg_test() {

  double tol = 1e-5;

  {
    SymGroup pg;
    Lattice::fcc().generate_point_group(pg, tol);
    BOOST_CHECK_EQUAL(pg.size(), 48);
  }
  {
    SymGroup pg;
    Lattice::bcc().generate_point_group(pg, tol);
    BOOST_CHECK_EQUAL(pg.size(), 48);
  }
  {
    SymGroup pg;
    Lattice::cubic().generate_point_group(pg, tol);
    BOOST_CHECK_EQUAL(pg.size(), 48);
  }
  {
    SymGroup pg;
    Lattice::hexagonal().generate_point_group(pg, tol);
    BOOST_CHECK_EQUAL(pg.size(), 24);
  }
}

void lattice_is_equivalent_test() {

  double tol = TOL;

  {
    SymGroup pg;
    Lattice fcc = Lattice::fcc();
    fcc.generate_point_group(pg, tol);

    for(const auto &op : pg) {
      BOOST_CHECK_EQUAL(fcc.is_equivalent(copy_apply(op, fcc), tol), 1);
    }
  }
  {
    SymGroup pg;
    Lattice bcc = Lattice::bcc();
    bcc.generate_point_group(pg, tol);

    for(const auto &op : pg) {
      BOOST_CHECK_EQUAL(bcc.is_equivalent(copy_apply(op, bcc), tol), 1);
    }
  }
  {
    SymGroup pg;
    Lattice cubic = Lattice::cubic();
    cubic.generate_point_group(pg, tol);

    for(const auto &op : pg) {
      BOOST_CHECK_EQUAL(cubic.is_equivalent(copy_apply(op, cubic), tol), 1);
    }
  }
  {
    SymGroup pg;
    Lattice hex = Lattice::hexagonal();
    hex.generate_point_group(pg, tol);

    for(const auto &op : pg) {
      BOOST_CHECK_EQUAL(hex.is_equivalent(copy_apply(op, hex), tol), 1);
    }
  }
}

void lattice_read_test() {
  std::stringstream lat_stream("   2.5000\n"
                               "   1.1 1.2 1.3\n"
                               "   1.4 1.5 1.6\n"
                               "   1.7 1.8 1.9\n");
  Lattice testlat;
  testlat.read(lat_stream);
  Eigen::Matrix3d latmat;
  latmat <<
         2.75, 3.5, 4.25,
               3.0, 3.75, 4.5,
               3.25, 4.0, 4.75;
  BOOST_CHECK(almost_equal(testlat.lat_column_mat(), latmat, 1e-8));

}

void lattice_superduper_test() {
  //This point group will remain empty so that it checks more cases
  SymGroup pg;
  Lattice lat(Lattice::fcc());

  ScelEnumProps enum_props(1, 6);
  SupercellEnumerator<Lattice> enumerator(lat, pg, enum_props);

  std::vector<Lattice> lat_list(enumerator.begin(), enumerator.end());

  for(auto it1 = lat_list.cbegin(); it1 != lat_list.cend(); ++it1) {
    for(auto it2 = it1 + 1; it2 != lat_list.cend(); ++it2) {
      Lattice sdlat = superdupercell(*it1, *it2);
      BOOST_CHECK(sdlat.is_supercell_of(*it1));
      BOOST_CHECK(sdlat.is_supercell_of(*it2));
    }
  }

}


BOOST_AUTO_TEST_SUITE(LatticeTest)

BOOST_AUTO_TEST_CASE(ReadTest) {
  lattice_read_test();
}

BOOST_AUTO_TEST_CASE(PointGroupTest) {
  lattice_pg_test();
}

BOOST_AUTO_TEST_CASE(IsEquivalentTest) {
  lattice_is_equivalent_test();
}

BOOST_AUTO_TEST_CASE(SuperDuperTest) {
  lattice_superduper_test();

}


BOOST_AUTO_TEST_SUITE_END()
