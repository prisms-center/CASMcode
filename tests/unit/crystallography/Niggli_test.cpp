#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/crystallography/Niggli.hh"

/// What is being used to test it:
#include "casm/container/LinearAlgebra.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SupercellEnumerator.hh"

namespace CASM {

  void confirm_lattice(const Lattice &known_niggli_form, const Eigen::Matrix3i &skewed_unimodular) {
    assert(skewed_unimodular.determinant() == 1);
    assert(is_niggli(known_niggli_form, CASM::TOL));

    //Lattice non_niggli(known_niggli_form.lat_column_mat()*skewed_unimodular.cast<double>());
    Lattice non_niggli(known_niggli_form.lat_column_mat()*skewed_unimodular.cast<double>());

    BOOST_CHECK(!is_niggli(non_niggli, CASM::TOL));

    Lattice reniggli = niggli(non_niggli, CASM::TOL);

    BOOST_CHECK(known_niggli_form == reniggli);

    BOOST_CHECK(niggli(niggli(non_niggli, CASM::TOL), CASM::TOL) == niggli(non_niggli, CASM::TOL));

    return;
  }


  void confirm_fcc_lattice(const Eigen::Matrix3i &skewed_unimodular) {
    BOOST_TEST_MESSAGE("Checking fcc lattice");
    CASM::confirm_lattice(CASM::Lattice::fcc(), skewed_unimodular);
    CASM::confirm_lattice(CASM::Lattice::fcc(), skewed_unimodular.transpose());
    return;
  }

  void confirm_bcc_lattice(const Eigen::Matrix3i &skewed_unimodular) {
    BOOST_TEST_MESSAGE("Checking bcc lattice");
    CASM::confirm_lattice(CASM::Lattice::bcc(), skewed_unimodular);
    CASM::confirm_lattice(CASM::Lattice::bcc(), skewed_unimodular.transpose());
    return;
  }

  void confirm_hexagonal_lattice(const Eigen::Matrix3i &skewed_unimodular) {
    BOOST_TEST_MESSAGE("Checking hexagonal lattice");
    CASM::confirm_lattice(CASM::Lattice::hexagonal(), skewed_unimodular);
    CASM::confirm_lattice(CASM::Lattice::hexagonal(), skewed_unimodular.transpose());
    return;
  }

  void confirm_cubic_lattice(const Eigen::Matrix3i &skewed_unimodular) {
    BOOST_TEST_MESSAGE("Checking cubic lattice");
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

    BOOST_CHECK(is_symmetric(symmat));
    BOOST_CHECK(!is_persymmetric(symmat));

    persymmat << 4, 3, 2, 1,
              7, 6, 5, 2,
              9, 8, 6, 3,
              10, 9, 7, 4;

    BOOST_CHECK(!is_symmetric(persymmat));
    BOOST_CHECK(is_persymmetric(persymmat));
  }

  //See issue #153 on github: https://github.com/prisms-center/CASMcode-dev/issues/153
  void single_dimension_test() {
    Lattice testlat = Lattice::fcc();
    SymGroup pg;
    testlat.generate_point_group(pg);

    int dims = 1;
    int minvol = 1;
    int maxvol = 10;

    SupercellEnumerator<Lattice> latenumerator(testlat, pg, minvol, maxvol + 1, dims);
    std::vector<Lattice> enumerated_lat(latenumerator.begin(), latenumerator.end());


    int l = 1;
    for(auto it = enumerated_lat.begin(); it != enumerated_lat.end(); ++it) {
      Eigen::Matrix3i comp_transmat;
      comp_transmat << (l), 0, 0,
                    0, 1, 0,
                    0, 0, 1;

      Lattice comparelat = make_supercell(testlat, comp_transmat);

      Lattice nigglicompare = canonical_equivalent_lattice(comparelat, pg, TOL);
      Lattice nigglitest = canonical_equivalent_lattice(*it, pg, TOL);

      BOOST_CHECK(nigglicompare == nigglitest);
      l++;
    }
    return;
  }
}

BOOST_AUTO_TEST_SUITE(NiggliTest)

BOOST_AUTO_TEST_CASE(SymmetricTest) {
  CASM::symmetric_testing();
}

BOOST_AUTO_TEST_CASE(EeasyTests) {
  Eigen::Matrix3i skewed_unimodular;
  skewed_unimodular << 1, 2, 3,
                    0, 1, 4,
                    0, 0, 1;

  CASM::confirm_fcc_lattice(skewed_unimodular);
  CASM::confirm_bcc_lattice(skewed_unimodular);
  CASM::confirm_cubic_lattice(skewed_unimodular);
  CASM::confirm_hexagonal_lattice(skewed_unimodular.transpose());
}

BOOST_AUTO_TEST_CASE(EvilNiggliTest) {
  CASM::single_dimension_test();
}


BOOST_AUTO_TEST_SUITE_END()
