#include "gtest/gtest.h"

/// What is being tested:
#include "casm/crystallography/Niggli.hh"
#include "casm/crystallography/CanonicalForm.hh"

/// What is being used to test it:
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/ScelEnum.hh"
#include "crystallography/TestStructures.hh"

namespace CASM {

  void confirm_lattice(const Lattice &known_niggli_form, const Eigen::Matrix3i &skewed_unimodular) {
    assert(skewed_unimodular.determinant() == 1);
    assert(is_niggli(known_niggli_form, CASM::TOL));

    Lattice non_niggli(known_niggli_form.lat_column_mat()*skewed_unimodular.cast<double>());

    EXPECT_TRUE(!is_niggli(non_niggli, CASM::TOL));

    Lattice reniggli = niggli(non_niggli, CASM::TOL);

    EXPECT_TRUE(is_niggli(known_niggli_form, CASM::TOL));
    EXPECT_TRUE(is_niggli(reniggli, CASM::TOL));
    EXPECT_TRUE(known_niggli_form == reniggli);
    EXPECT_TRUE(niggli(niggli(non_niggli, CASM::TOL), CASM::TOL) == niggli(non_niggli, CASM::TOL));


    return;
  }


  void confirm_fcc_lattice(const Eigen::Matrix3i &skewed_unimodular) {
    /* BOOST_TEST_MESSAGE("Checking fcc lattice"); */
    CASM::confirm_lattice(CASM::Lattice::fcc(), skewed_unimodular);
    CASM::confirm_lattice(CASM::Lattice::fcc(), skewed_unimodular.transpose());
    return;
  }

  void confirm_bcc_lattice(const Eigen::Matrix3i &skewed_unimodular) {
    /* BOOST_TEST_MESSAGE("Checking bcc lattice"); */
    CASM::confirm_lattice(CASM::Lattice::bcc(), skewed_unimodular);
    CASM::confirm_lattice(CASM::Lattice::bcc(), skewed_unimodular.transpose());
    return;
  }

  void confirm_hexagonal_lattice(const Eigen::Matrix3i &skewed_unimodular) {
    /* BOOST_TEST_MESSAGE("Checking hexagonal lattice"); */
    CASM::confirm_lattice(CASM::Lattice::hexagonal(), skewed_unimodular);
    CASM::confirm_lattice(CASM::Lattice::hexagonal(), skewed_unimodular.transpose());
    return;
  }

  void confirm_cubic_lattice(const Eigen::Matrix3i &skewed_unimodular) {
    /* BOOST_TEST_MESSAGE("Checking cubic lattice"); */
    CASM::confirm_lattice(CASM::Lattice::cubic(), skewed_unimodular);
    CASM::confirm_lattice(CASM::Lattice::cubic(), skewed_unimodular.transpose());
    return;
  }

  //Test for checking whether matrices are symmetric or persymmetric
  void symmetric_testing() {
    Eigen::MatrixXd symmat(5, 5), persymmat(4, 4), nonsymmat;

    symmat << 1, 2, 3, 4, 5,
           2, 6, 7, 8, 9,
           3, 7, 10, 11, 12,
           4, 8, 11, 13, 14,
           5, 9, 12, 14, 15;

    EXPECT_TRUE(is_symmetric(symmat));
    EXPECT_TRUE(!is_persymmetric(symmat));

    persymmat << 4, 3, 2, 1,
              7, 6, 5, 2,
              9, 8, 6, 3,
              10, 9, 7, 4;

    EXPECT_TRUE(!is_symmetric(persymmat));
    EXPECT_TRUE(is_persymmetric(persymmat));
  }

  //See issue #153 on github: https://github.com/prisms-center/CASMcode-dev/issues/153
  void single_dimension_test() {
    Lattice testlat = Lattice::fcc();
    SymGroup pg = SymGroup::lattice_point_group(testlat);

    std::string dirs = "a";
    int minvol = 1;
    int maxvol = 10;

    ScelEnumProps enum_props(minvol, maxvol + 1, dirs);
    SuperlatticeEnumerator latenumerator(pg.begin(), pg.end(), testlat, enum_props);
    std::vector<Lattice> enumerated_lat(latenumerator.begin(), latenumerator.end());


    int l = 1;
    for(auto it = enumerated_lat.begin(); it != enumerated_lat.end(); ++it) {
      Eigen::Matrix3i comp_transmat;
      comp_transmat << (l), 0, 0,
                    0, 1, 0,
                    0, 0, 1;

      Lattice comparelat = make_superlattice(testlat, comp_transmat);

      Lattice nigglicompare = xtal::canonical::equivalent(comparelat, pg, TOL);
      Lattice nigglitest = xtal::canonical::equivalent(*it, pg, TOL);

      EXPECT_TRUE(nigglicompare == nigglitest);
      l++;
    }
    return;
  }

  void standard_orientation_compare_test() {

    /* BOOST_TEST_MESSAGE("Checking standard_orientation_compare"); */

    double tol = TOL;

    //This is a known supercell of ZrO

    Eigen::Matrix3d lat_mat_A;
    lat_mat_A << 3.233986860000,  0.000000000000,  0.000000000000,
              0.000000000000,  0.000000000000,  5.601429540000,
              0.000000000000, -5.168678340000,  0.000000000000;

    Lattice lat_A(lat_mat_A);

    Eigen::Matrix3d lat_mat_A2;
    lat_mat_A2 << 3.233986860000,  0.000000000000,  0.000000000000,
               2.22045e-16,     0.000000000000,  5.601429540000,
               0.000000000000, -5.168678340000,  0.000000000000;

    Lattice lat_A2(lat_mat_A2);

    Eigen::Matrix3d lat_mat_B;
    lat_mat_B << 3.233986860000,  0.000000000000,  0.000000000000,
              0.000000000000,  0.000000000000, -5.601429540000,
              0.000000000000,  5.168678340000,  0.000000000000;

    Lattice lat_B(lat_mat_B);

    EXPECT_EQ(xtal::standard_orientation_compare(lat_mat_A, lat_mat_B, tol), true);
    EXPECT_EQ(xtal::standard_orientation_compare(lat_mat_B, lat_mat_A, tol), false);

    EXPECT_EQ(xtal::standard_orientation_compare(lat_mat_A2, lat_mat_B, tol), true);
    EXPECT_EQ(xtal::standard_orientation_compare(lat_mat_B, lat_mat_A2, tol), false);

    EXPECT_EQ(xtal::standard_orientation_compare(lat_mat_A, lat_mat_A2, tol), false);
    EXPECT_EQ(xtal::standard_orientation_compare(lat_mat_A2, lat_mat_A, tol), false);

    Structure prim(test::ZrO_prim());
    Lattice canon_A = xtal::canonical::equivalent(lat_A, prim.point_group(), tol);
    Lattice canon_A2 = xtal::canonical::equivalent(lat_A2, prim.point_group(), tol);
    Lattice canon_B = xtal::canonical::equivalent(lat_B, prim.point_group(), tol);

    EXPECT_EQ(canon_A == canon_A2, true);
    EXPECT_EQ(canon_A2 == canon_B, true);
    EXPECT_EQ(canon_A == canon_B, true);

  }
}

TEST(NiggliTest, SymmetricTest) {
  CASM::symmetric_testing();
}

TEST(NiggliTest, EasyTests) {
  Eigen::Matrix3i skewed_unimodular;
  skewed_unimodular << 1, 2, 3,
                    0, 1, 4,
                    0, 0, 1;

  CASM::confirm_fcc_lattice(skewed_unimodular);
  CASM::confirm_bcc_lattice(skewed_unimodular);
  CASM::confirm_cubic_lattice(skewed_unimodular);
  CASM::confirm_hexagonal_lattice(skewed_unimodular);
}

TEST(NiggliTest, EvilNiggliTest) {
  CASM::single_dimension_test();
}

