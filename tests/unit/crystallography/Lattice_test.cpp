#include "gtest/gtest.h"

/// What is being tested:
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/LatticeIsEquivalent.hh"

/// What is being used to test it:
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/symmetry/SymGroup.hh"

using namespace CASM;
using xtal::SuperlatticeEnumerator;
using xtal::ScelEnumProps;

void lattice_pg_test() {

  double tol = 1e-5;

  {
    EXPECT_EQ(SymGroup::lattice_point_group(Lattice::fcc()).size(), 48);

    EXPECT_EQ(SymGroup::lattice_point_group(Lattice::bcc()).size(), 48);

    EXPECT_EQ(SymGroup::lattice_point_group(Lattice::cubic()).size(), 48);

    EXPECT_EQ(SymGroup::lattice_point_group(Lattice::hexagonal()).size(), 24);
  }
}

void lattice_is_equivalent_test() {

  {
    Lattice fcc = Lattice::fcc();
    SymGroup pg = SymGroup::lattice_point_group(fcc);
    for(const auto &op : pg) {
      EXPECT_TRUE(xtal::is_equivalent(fcc, sym::copy_apply(op, fcc)));
    }
  }
  {

    Lattice bcc = Lattice::bcc();
    SymGroup pg = SymGroup::lattice_point_group(bcc);

    for(const auto &op : pg) {
      EXPECT_TRUE(xtal::is_equivalent(bcc, sym::copy_apply(op, bcc)));
    }
  }
  {

    Lattice cubic = Lattice::cubic();
    SymGroup pg = SymGroup::lattice_point_group(cubic);

    for(const auto &op : pg) {
      EXPECT_TRUE(xtal::is_equivalent(cubic, sym::copy_apply(op, cubic)));
    }
  }
  {
    Lattice hex = Lattice::hexagonal();
    SymGroup pg = SymGroup::lattice_point_group(hex);

    for(const auto &op : pg) {
      EXPECT_TRUE(xtal::is_equivalent(hex, sym::copy_apply(op, hex)));
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
  EXPECT_TRUE(almost_equal(testlat.lat_column_mat(), latmat, 1e-8));

}

void lattice_superduper_test() {
  //This point group will remain empty so that it checks more cases
  SymGroup pg;
  Lattice lat(Lattice::fcc());

  ScelEnumProps enum_props(1, 6);
  SuperlatticeEnumerator enumerator(pg.begin(), pg.end(), lat, enum_props);

  std::vector<Lattice> lat_list(enumerator.begin(), enumerator.end());

  for(auto it1 = lat_list.cbegin(); it1 != lat_list.cend(); ++it1) {
    for(auto it2 = it1 + 1; it2 != lat_list.cend(); ++it2) {
      Lattice sdlat = xtal::make_superduperlattice(*it1, *it2);
      EXPECT_TRUE(xtal::is_superlattice(sdlat, *it1, sdlat.tol()).first);
      EXPECT_TRUE(xtal::is_superlattice(sdlat, *it2, sdlat.tol()).first);
    }
  }

}



TEST(LatticeTest, ReadTest) {
  lattice_read_test();
}

TEST(LatticeTest, PointGroupTest) {
  lattice_pg_test();
}

TEST(LatticeTest, IsEquivalentTest) {
  lattice_is_equivalent_test();
}

TEST(LatticeTest, SuperDuperTest) {
  lattice_superduper_test();

}
