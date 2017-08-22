#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/Configuration.hh"

/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/symmetry/SymInfo.hh"
#include "casm/app/AppIO.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(ConfigurationTest)

BOOST_AUTO_TEST_CASE(Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.get_prim().lattice().vectors();

  Supercell scel {&primclex, Lattice(a, b, c)};

  Configuration config(scel);
  BOOST_CHECK_EQUAL(config.size(), 1);

  // set occupation
  BOOST_CHECK_EQUAL(config.has_occupation(), false);

  config.init_occupation();
  BOOST_CHECK_EQUAL(config.has_occupation(), true);

  config.clear_occupation();
  BOOST_CHECK_EQUAL(config.has_occupation(), false);

  config.set_occupation(Array<int>(1, 0));
  BOOST_CHECK_EQUAL(config.has_occupation(), true);

  for(int i = 0; i < 3; ++i) {
    config.set_occ(0, i);
    BOOST_CHECK_EQUAL(config.occ(0), i);
  }

  // set displacement
  typedef Configuration::displacement_matrix_t disp_matrix_t;
  BOOST_CHECK_EQUAL(config.has_displacement(), false);

  config.init_displacement();
  BOOST_CHECK_EQUAL(config.has_displacement(), true);

  config.clear_displacement();
  BOOST_CHECK_EQUAL(config.has_displacement(), false);

  config.set_displacement(disp_matrix_t::Zero(3, 1));
  BOOST_CHECK_EQUAL(config.has_displacement(), true);

  for(int i = 0; i < 3; ++i) {
    Eigen::Vector3d d = Eigen::Vector3d::Zero();
    d(i) = 0.001;
    config.set_disp(0, d);
    BOOST_CHECK_EQUAL(almost_equal(config.disp(0), d), true);
  }

  // set deformation
  BOOST_CHECK_EQUAL(config.has_deformation(), false);

  config.init_deformation();
  BOOST_CHECK_EQUAL(config.has_deformation(), true);

  config.clear_deformation();
  BOOST_CHECK_EQUAL(config.has_deformation(), false);

  config.set_deformation(Eigen::Matrix3d::Zero());
  BOOST_CHECK_EQUAL(config.has_deformation(), true);

}

BOOST_AUTO_TEST_CASE(Test2) {

  // test Configuration::fill_supercell

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.get_prim().lattice().vectors();

  Supercell scel {&primclex, Lattice(c, a - b, a + b - c)};

  Configuration config(scel);
  BOOST_CHECK_EQUAL(config.size(), 2);

  // include occupation only
  config.set_occupation({1, 0});

  {
    // Identity op
    Supercell scel {&primclex, Lattice(c, a - b, a + b - c)};
    const SymGroup &fg = config.get_supercell().factor_group();
    Configuration filled = config.fill_supercell(scel, fg[0]);

    Configuration check(scel);
    check.set_occupation({1, 0});

    BOOST_CHECK_EQUAL(filled, check);
  }

  {
    // supercell
    Supercell scel {&primclex, Lattice(c, a - b, 2.*(a + b - c))};
    const SymGroup &fg = config.get_supercell().factor_group();
    Configuration filled = config.fill_supercell(scel, fg[0]);

    Configuration check(scel);
    check.set_occupation({1, 0, 1, 0});

    BOOST_CHECK_EQUAL(filled, check);
  }

  {
    // 90 deg rotation
    Supercell scel {&primclex, Lattice(c, a - b, a + b - c)};
    const SymGroup &fg = config.get_supercell().factor_group();
    //std::cout << to_string(fg.info(1), CART) << std::endl;
    Configuration filled = config.fill_supercell(scel, fg[1]);

    Configuration check(scel);
    check.set_occupation({1, 0});

    BOOST_CHECK_EQUAL(filled, check);
  }


  // include occupation & displacements
  Eigen::Vector3d dzero(0., 0., 0.);
  Eigen::Vector3d dx(0.001, 0., 0.);
  Eigen::Vector3d dy(0., 0.001, 0.);
  Eigen::Vector3d dz(0., 0., 0.001);

  config.init_displacement();
  config.set_disp(0, dx);

  {
    // Identity op
    Supercell scel {&primclex, Lattice(c, a - b, a + b - c)};
    const SymGroup &fg = config.get_supercell().factor_group();
    Configuration filled = config.fill_supercell(scel, fg[0]);

    Configuration check(scel);
    check.set_occupation({1, 0});
    check.init_displacement();
    check.set_disp(0, dx);


    BOOST_CHECK_EQUAL(filled, check);
  }

  {
    // supercell
    Supercell scel {&primclex, Lattice(c, a - b, 2.*(a + b - c))};
    const SymGroup &fg = config.get_supercell().factor_group();
    Configuration filled = config.fill_supercell(scel, fg[0]);

    Configuration check(scel);
    check.set_occupation({1, 0, 1, 0});
    check.init_displacement();
    check.set_disp(0, dx);
    check.set_disp(2, dx);

    BOOST_CHECK_EQUAL(filled, check);
  }

  {
    // 90 deg rotation
    Supercell scel {&primclex, Lattice(c, a - b, a + b - c)};
    const SymGroup &fg = config.get_supercell().factor_group();
    Configuration filled = config.fill_supercell(scel, fg[1]);

    Configuration check(scel);
    check.set_occupation({1, 0});
    check.init_displacement();
    check.set_disp(0, dy);

    BOOST_CHECK_EQUAL(filled, check);
  }

}

BOOST_AUTO_TEST_SUITE_END()
