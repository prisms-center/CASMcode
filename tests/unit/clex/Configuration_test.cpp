#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/Configuration_impl.hh"

/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/app/AppIO.hh"
#include "casm/crystallography/Structure.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(ConfigurationTest)

BOOST_AUTO_TEST_CASE(Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  Supercell scel {&primclex, Lattice(a, b, c)};

  Configuration config(scel);
  BOOST_CHECK_EQUAL(config.size(), 1);

  // set occupation
  BOOST_CHECK_EQUAL(config.has_occupation(), false);

  config.init_occupation();
  BOOST_CHECK_EQUAL(config.has_occupation(), true);

  config.clear_occupation();
  BOOST_CHECK_EQUAL(config.has_occupation(), false);

  config.set_occupation(std::vector<int>({0}));
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
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  Supercell scel {&primclex, Lattice(c, a - b, a + b - c)};

  Configuration config(scel);
  BOOST_CHECK_EQUAL(config.size(), 2);

  // include occupation only
  config.set_occupation({1, 0});

  {
    // Identity op
    Supercell scel {&primclex, Lattice(c, a - b, a + b - c)};
    const SymGroup &fg = config.supercell().factor_group();
    Configuration filled = config.fill_supercell(scel, fg[0]);

    Configuration check(scel);
    check.set_occupation({1, 0});

    BOOST_CHECK_EQUAL(filled, check);
  }

  {
    // supercell
    Supercell scel {&primclex, Lattice(c, a - b, 2.*(a + b - c))};
    const SymGroup &fg = config.supercell().factor_group();
    Configuration filled = config.fill_supercell(scel, fg[0]);

    Configuration check(scel);
    check.set_occupation({1, 0, 1, 0});

    BOOST_CHECK_EQUAL(filled, check);
  }

  {
    // 90 deg rotation
    Supercell scel {&primclex, Lattice(c, a - b, a + b - c)};
    const SymGroup &fg = config.supercell().factor_group();
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
    const SymGroup &fg = config.supercell().factor_group();
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
    const SymGroup &fg = config.supercell().factor_group();
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
    const SymGroup &fg = config.supercell().factor_group();
    Configuration filled = config.fill_supercell(scel, fg[1]);

    Configuration check(scel);
    check.set_occupation({1, 0});
    check.init_displacement();
    check.set_disp(0, dy);

    BOOST_CHECK_EQUAL(filled, check);
  }

}

BOOST_AUTO_TEST_CASE(Test3) {
  // test ConfigCanonicalForm functions

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  {
    // supercell (standard cubic FCC)
    Supercell scel {&primclex, Lattice(a - b + c, a + b - c, b + c - a)};
    std::cout << scel.lattice().lat_column_mat() << std::endl;
    BOOST_CHECK_EQUAL(scel.is_canonical(), true);

    {
      Configuration config(scel);
      config.set_occupation({1, 0, 0, 0});
      BOOST_CHECK_EQUAL(config.is_canonical(), true);
      BOOST_CHECK_EQUAL(config.is_primitive(), true);
      BOOST_CHECK_EQUAL(config.invariant_subgroup().size(), 48);

      {
        Configuration test(scel);
        test.set_occupation({1, 0, 0, 0});
        BOOST_CHECK_EQUAL(config == test, true);
        BOOST_CHECK_EQUAL(config.is_equivalent(test), true);
        BOOST_CHECK_EQUAL(test.is_equivalent(config), true);
        BOOST_CHECK_EQUAL(test < config, false);
        BOOST_CHECK_EQUAL(config < test, false);
      }

      {
        Configuration test(scel);
        test.set_occupation({0, 1, 0, 0});
        BOOST_CHECK_EQUAL(config == test, false);
        BOOST_CHECK_EQUAL(config.is_equivalent(test), true);
        BOOST_CHECK_EQUAL(test.is_equivalent(config), true);
        BOOST_CHECK_EQUAL(test < config, true);
        BOOST_CHECK_EQUAL(config < test, false);

        auto to_canonical = test.to_canonical();
        BOOST_CHECK_EQUAL(to_canonical.factor_group_index(), 0);
        BOOST_CHECK_EQUAL(to_canonical.translation_index(), 0);
        BOOST_CHECK_EQUAL(copy_apply(to_canonical, test) == config, true);
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
