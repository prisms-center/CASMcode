#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/Configuration_impl.hh"

/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/app/AppIO.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/database/ConfigDatabase.hh"

using namespace CASM;

/// Check that name -> configuration recreates an equivalent configuration
void check_made_from_name(const Configuration &config, std::string name) {
  BOOST_CHECK_EQUAL(true, true);

  // make configuration from name and database
  Configuration made_from_name = make_configuration(config.primclex(), name);
  BOOST_CHECK_EQUAL(true, true);

  // check that Supercell are equivalent (though vectors may be different)
  BOOST_CHECK_EQUAL(
    config.supercell().lattice().is_equivalent(made_from_name.supercell().lattice()),
    true);

  // fill supercell and check that config are identical
  Configuration in_same_supercell = made_from_name.fill_supercell(config.supercell());
  BOOST_CHECK_EQUAL((config == in_same_supercell), true);
}

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

BOOST_AUTO_TEST_CASE(TestConfigurationName) {
  // test Configuration::generate_name_impl
  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  auto &db = primclex.db<Configuration>();

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  {
    // prim cell
    Supercell scel {&primclex, Lattice(a, b, c)};

    {
      // canonical scel, canonical primitive occ
      Configuration config(scel);
      config.set_occupation({0});

      // not in datbase -> id == "none"
      BOOST_CHECK_EQUAL(config.name(), "SCEL1_1_1_1_0_0_0/none");

      auto res = db.insert(config);

      // not from database, but now in it -> id == "0"
      BOOST_CHECK_EQUAL(config.name(), "SCEL1_1_1_1_0_0_0/0");

      // in datbase -> id == "0"
      BOOST_CHECK_EQUAL(res.first->name(), "SCEL1_1_1_1_0_0_0/0");

      check_made_from_name(config, res.first->name());
    }

    {
      // canonical scel, canonical primitive occ
      Configuration config(scel);
      config.set_occupation({1});
      // not in datbase -> id == "none"
      BOOST_CHECK_EQUAL(config.name(), "SCEL1_1_1_1_0_0_0/none");

      auto res = db.insert(config);

      // not from database, but now in it -> id == "1"
      BOOST_CHECK_EQUAL(config.name(), "SCEL1_1_1_1_0_0_0/1");

      // in datbase -> id == "1"
      BOOST_CHECK_EQUAL(res.first->name(), "SCEL1_1_1_1_0_0_0/1");

      check_made_from_name(config, res.first->name());
    }

    {
      // canonical scel, canonical primitive occ
      Configuration config(scel);
      config.set_occupation({0});

      // not from database, but now in it -> id == "0"
      BOOST_CHECK_EQUAL(config.name(), "SCEL1_1_1_1_0_0_0/0");
    }

    while(db.size()) {
      db.erase(db.begin());
    }
  }

  {
    // standard cubic FCC unit cell
    Supercell scel {&primclex, Lattice(c + b - a, a - b + c, a + b - c)};

    {
      // canonical scel, canonical primitive occ
      Configuration config(scel);
      config.set_occupation({1, 0, 0, 0});
      BOOST_CHECK_EQUAL(config.name(), "SCEL4_2_2_1_1_1_0/none");

      auto res = db.insert(config.in_canonical_supercell());

      BOOST_CHECK_EQUAL(res.first->name(), "SCEL4_2_2_1_1_1_0/0");

      check_made_from_name(config, "SCEL4_2_2_1_1_1_0/0");
    }

    {
      // canonical scel, non-canonical primitive occ
      Configuration config(scel);
      config.set_occupation({0, 1, 0, 0});

      BOOST_CHECK_EQUAL(config.name(), "SCEL4_2_2_1_1_1_0/0.equiv.0.1");

      auto res = db.insert(config.in_canonical_supercell());
      BOOST_CHECK_EQUAL(res.second, false);

      BOOST_CHECK_EQUAL(config.name(), "SCEL4_2_2_1_1_1_0/0.equiv.0.1");

      BOOST_CHECK_EQUAL(res.first->name(), "SCEL4_2_2_1_1_1_0/0");

      check_made_from_name(config, "SCEL4_2_2_1_1_1_0/0.equiv.0.1");
    }

    {
      // canonical scel, canonical non-primitive occ
      Configuration config(scel);
      config.set_occupation({0, 0, 0, 0});

      BOOST_CHECK_EQUAL(config.name(), "SCEL4_2_2_1_1_1_0/none");

      auto res = db.insert(config.in_canonical_supercell());
      BOOST_CHECK_EQUAL(res.second, true);

      // primitive does not yet exist in database
      BOOST_CHECK_EQUAL(config.name(), "SCEL4_2_2_1_1_1_0/1");
      BOOST_CHECK_EQUAL(res.first->name(), "SCEL4_2_2_1_1_1_0/1");

      // insert primitive
      res = db.insert(config.primitive().in_canonical_supercell());

      // primitive does exist in database
      // - it gets index 2, even though it is the only config from this supercell
      //   in the database, because we don't re-use indices
      BOOST_CHECK_EQUAL(res.first->name(), "SCEL1_1_1_1_0_0_0/2");

      // make config from primitive name
      check_made_from_name(config, "SCEL4_2_2_1_1_1_0/super.0.SCEL1_1_1_1_0_0_0/2.equiv.0.0");
    }
  }

  {
    // non-canonical FCC unit cell, equivalent to standard FCC unit cell
    Eigen::Vector3d standard_a = c + b - a;
    Eigen::Vector3d standard_b = a - b + c;
    Eigen::Vector3d standard_c = a + b - c;

    Supercell scel {&primclex, Lattice(standard_a, standard_b, standard_c + standard_b)};

    {
      // non-canonical scel, canonical primitive occ
      Configuration config(scel);
      config.set_occupation({1, 0, 0, 0});

      // having a different, but equivalent supercell, should not change name ^ see above
      BOOST_CHECK_EQUAL(config.name(), "SCEL4_2_2_1_1_1_0/0");

      auto res = db.insert(config.in_canonical_supercell());
      BOOST_CHECK_EQUAL(res.second, false);
      BOOST_CHECK_EQUAL(config.name(), "SCEL4_2_2_1_1_1_0/0");
      BOOST_CHECK_EQUAL(res.first->name(), "SCEL4_2_2_1_1_1_0/0");

      check_made_from_name(config, "SCEL4_2_2_1_1_1_0/0");
    }

    {
      // non-canonical scel, non-canonical primitive occ
      Configuration config(scel);
      config.set_occupation({0, 1, 0, 0});

      // having a different, but equivalent supercell, should not change name ^ see above
      BOOST_CHECK_EQUAL(config.name(), "SCEL4_2_2_1_1_1_0/0.equiv.0.3");

      auto res = db.insert(config.in_canonical_supercell());
      BOOST_CHECK_EQUAL(res.second, false);

      BOOST_CHECK_EQUAL(config.name(), "SCEL4_2_2_1_1_1_0/0.equiv.0.3");

      BOOST_CHECK_EQUAL(res.first->name(), "SCEL4_2_2_1_1_1_0/0");

      check_made_from_name(config, "SCEL4_2_2_1_1_1_0/0.equiv.0.3");
    }

    {
      // equivalent, but non-canonical scel, canonical non-primitive occ
      Configuration config(scel);
      config.set_occupation({0, 0, 0, 0});

      // having a different, but equivalent supercell, should not change name ^ see above
      BOOST_CHECK_EQUAL(config.name(), "SCEL4_2_2_1_1_1_0/1");

      auto res = db.insert(config.in_canonical_supercell());
      BOOST_CHECK_EQUAL(res.second, false);

      check_made_from_name(config, "SCEL4_2_2_1_1_1_0/super.0.SCEL1_1_1_1_0_0_0/2.equiv.0.0");

    }
  }

  {
    // Test non-canonical equivalent unit cell
    Eigen::Vector3d standard_a = c + b - a;
    Eigen::Vector3d standard_b = a - b + c;
    Eigen::Vector3d standard_c = a + b - c;

    Lattice canon_lat = canonical_equivalent_lattice(
                          Lattice(standard_a, standard_b, 2 * standard_c),
                          primclex.prim().point_group(),
                          primclex.crystallography_tol());
    Lattice test_lat;
    Index scel_op_index = 0;
    for(const auto &op : primclex.prim().point_group()) {
      test_lat = copy_apply(op, canon_lat);
      if(!test_lat.is_equivalent(canon_lat)) {
        break;
      }
      ++scel_op_index;
    }

    Supercell scel {&primclex, test_lat};

    {
      // non-canonical, non-equivalent, scel, canonical primitive occ
      Configuration config(scel);
      config.set_occupation({1, 0, 0, 0, 0, 0, 0, 0});

      // having a different, but equivalent supercell, should not change name ^ see above
      BOOST_CHECK_EQUAL(config.name(), "SCEL8_4_2_1_1_3_2.4/super.1.SCEL8_4_2_1_1_3_2/none.equiv.0.0");

      auto res = db.insert(config.in_canonical_supercell());
      BOOST_CHECK_EQUAL(res.second, true);
      BOOST_CHECK_EQUAL(config.name(), "SCEL8_4_2_1_1_3_2.4/super.1.SCEL8_4_2_1_1_3_2/0.equiv.0.0");
      BOOST_CHECK_EQUAL(res.first->name(), "SCEL8_4_2_1_1_3_2/0");

      check_made_from_name(config, "SCEL8_4_2_1_1_3_2.4/super.1.SCEL8_4_2_1_1_3_2/0.equiv.0.0");
    }

    {
      // non-canonical, non-equivalent, scel, non-canonical primitive occ
      Configuration config(scel);
      config.set_occupation({0, 1, 0, 0, 0, 0, 0, 0});

      // having a different, but equivalent supercell, should not change name ^ see above
      BOOST_CHECK_EQUAL(config.name(), "SCEL8_4_2_1_1_3_2.4/super.1.SCEL8_4_2_1_1_3_2/0.equiv.0.1");

      auto res = db.insert(config.in_canonical_supercell());
      BOOST_CHECK_EQUAL(res.second, false);

      BOOST_CHECK_EQUAL(config.name(), "SCEL8_4_2_1_1_3_2.4/super.1.SCEL8_4_2_1_1_3_2/0.equiv.0.1");
      BOOST_CHECK_EQUAL(res.first->name(), "SCEL8_4_2_1_1_3_2/0");

      check_made_from_name(config, "SCEL8_4_2_1_1_3_2.4/super.1.SCEL8_4_2_1_1_3_2/0.equiv.0.1");
    }

    {
      // non-canonical, non-equivalent, scel, canonical non-primitive occ
      Configuration config(scel);
      config.set_occupation({0, 0, 0, 0, 0, 0, 0, 0});

      // having a different, but equivalent supercell, should not change name ^ see above
      BOOST_CHECK_EQUAL(config.name(), "SCEL8_4_2_1_1_3_2.4/super.0.SCEL1_1_1_1_0_0_0/2.equiv.0.0");

      auto res = db.insert(config.in_canonical_supercell());
      BOOST_CHECK_EQUAL(res.second, true);

      check_made_from_name(config, "SCEL8_4_2_1_1_3_2.4/super.0.SCEL1_1_1_1_0_0_0/2.equiv.0.0");
    }
  }

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
  // a = 0, 2, 2; b = 2, 0, 2; c = 2, 2, 0


  {
    // supercell (standard cubic FCC)
    Supercell scel {&primclex, Lattice(b + c - a, a + c - b, a + b - c)};
    //std::cout << scel.lattice().lat_column_mat() << std::endl;
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
        BOOST_CHECK_EQUAL(config.is_sym_equivalent(test), true);
        BOOST_CHECK_EQUAL(test.is_sym_equivalent(config), true);
        BOOST_CHECK_EQUAL(test < config, false);
        BOOST_CHECK_EQUAL(config < test, false);
      }

      {
        Configuration test(scel);
        test.set_occupation({0, 1, 0, 0});
        BOOST_CHECK_EQUAL(config == test, false);
        BOOST_CHECK_EQUAL(config.is_sym_equivalent(test), true);
        BOOST_CHECK_EQUAL(test.is_sym_equivalent(config), true);
        BOOST_CHECK_EQUAL(test < config, true);
        BOOST_CHECK_EQUAL(config < test, false);

        auto to_canonical = test.to_canonical();
        BOOST_CHECK_EQUAL(to_canonical.factor_group_index(), 0);
        BOOST_CHECK_EQUAL(to_canonical.translation_index(), 1);
        BOOST_CHECK_EQUAL(copy_apply(to_canonical, test) == config, true);
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
