#include "gtest/gtest.h"

/// What is being tested:
#include "casm/clex/Configuration_impl.hh"

/// What is being used to test it:

#include "Common.hh"
#include "TestConfiguration.hh"
#include "FCCTernaryProj.hh"
#include "casm/app/AppIO.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/database/ConfigDatabase.hh"

using namespace CASM;

/// Check that name -> configuration recreates an equivalent configuration
void check_made_from_name(const Configuration &config, std::string name) {
  EXPECT_EQ(true, true);

  // make configuration from name and database
  Configuration made_from_name = make_configuration(config.primclex(), name);
  EXPECT_EQ(true, true);

  // check that Supercell are equivalent (though vectors may be different)
  EXPECT_EQ(
    xtal::is_equivalent(config.supercell().lattice(), made_from_name.supercell().lattice()),
    true);

  // fill supercell and check that config are identical
  Configuration in_same_supercell = made_from_name.fill_supercell(config.supercell());
  EXPECT_EQ((config == in_same_supercell), true);
}

Lattice non_canonical_equiv_test_lat(const PrimClex &primclex) {

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  // Test non-canonical equivalent unit cell
  Eigen::Vector3d standard_a = c + b - a;
  Eigen::Vector3d standard_b = a - b + c;
  Eigen::Vector3d standard_c = a + b - c;

  Lattice canon_lat = xtal::canonical::equivalent(
                        Lattice(standard_a, standard_b, 2 * standard_c),
                        primclex.prim().point_group(),
                        primclex.crystallography_tol());
  Lattice test_lat;
  Index scel_op_index = 0;
  for(const auto &op : primclex.prim().point_group()) {
    test_lat = copy_apply(op, canon_lat);
    if(!xtal::is_equivalent(test_lat, canon_lat)) {
      return test_lat;
    }
    ++scel_op_index;
  }

  throw std::runtime_error("could not find non_canonical_equiv_test_lat");
}

TEST(ConfigurationTest, Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  Supercell scel {&primclex, Lattice(a, b, c)};

  Configuration config(scel);
  EXPECT_EQ(config.size(), 1);

  // set occupation
  EXPECT_EQ(config.occupation().size(), config.size());

  config.init_occupation();
  EXPECT_EQ(config.has_occupation(), true);

  config.set_occupation(std::vector<int>({0}));
  EXPECT_EQ(config.has_occupation(), true);

  for(int i = 0; i < 3; ++i) {
    config.set_occ(0, i);
    EXPECT_EQ(config.occ(0), i);
  }

  // set displacement
  /*
  typedef Configuration::displacement_matrix_t disp_matrix_t;
  EXPECT_EQ(config.has_displacement(), false);

  config.init_displacement();
  EXPECT_EQ(config.has_displacement(), true);

  config.clear_displacement();
  EXPECT_EQ(config.has_displacement(), false);

  config.set_displacement(disp_matrix_t::Zero(3, 1));
  EXPECT_EQ(config.has_displacement(), true);

  for(int i = 0; i < 3; ++i) {
    Eigen::Vector3d d = Eigen::Vector3d::Zero();
    d(i) = 0.001;
    config.set_disp(0, d);
    EXPECT_EQ(almost_equal(config.disp(0), d), true);
  }

  // set deformation
  EXPECT_EQ(config.has_deformation(), false);

  config.init_deformation();
  EXPECT_EQ(config.has_deformation(), true);

  config.clear_deformation();
  EXPECT_EQ(config.has_deformation(), false);

  config.set_deformation(Eigen::Matrix3d::Zero());
  EXPECT_EQ(config.has_deformation(), true);
  */
}

TEST(ConfigurationTest, TestConfigurationName) {
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
      config.set_occupation(std::vector<int>({0}));

      // not in datbase -> id == "none"
      EXPECT_EQ(config.name(), "SCEL1_1_1_1_0_0_0/none");

      auto res = db.insert(config);

      // not from database, but now in it -> id == "0"
      EXPECT_EQ(config.name(), "SCEL1_1_1_1_0_0_0/0");

      // in datbase -> id == "0"
      EXPECT_EQ(res.first->name(), "SCEL1_1_1_1_0_0_0/0");

      check_made_from_name(config, res.first->name());
    }

    {
      // canonical scel, canonical primitive occ
      Configuration config(scel);
      config.set_occupation(std::vector<int>({1}));
      // not in datbase -> id == "none"
      EXPECT_EQ(config.name(), "SCEL1_1_1_1_0_0_0/none");

      auto res = db.insert(config);

      // not from database, but now in it -> id == "1"
      EXPECT_EQ(config.name(), "SCEL1_1_1_1_0_0_0/1");

      // in datbase -> id == "1"
      EXPECT_EQ(res.first->name(), "SCEL1_1_1_1_0_0_0/1");

      check_made_from_name(config, res.first->name());
    }

    {
      // canonical scel, canonical primitive occ
      Configuration config(scel);
      config.set_occupation(std::vector<int>({0}));

      // not from database, but now in it -> id == "0"
      EXPECT_EQ(config.name(), "SCEL1_1_1_1_0_0_0/0");
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
      config.set_occupation(std::vector<int>({1, 0, 0, 0}));
      EXPECT_EQ(config.name(), "SCEL4_2_2_1_1_1_0/none");

      auto res = db.insert(config.in_canonical_supercell());

      EXPECT_EQ(res.first->name(), "SCEL4_2_2_1_1_1_0/0");

      check_made_from_name(config, "SCEL4_2_2_1_1_1_0/0");
    }

    {
      // canonical scel, non-canonical primitive occ
      Configuration config(scel);
      config.set_occupation(std::vector<int>({0, 1, 0, 0}));

      EXPECT_EQ(config.name(), "SCEL4_2_2_1_1_1_0/0.equiv.0.1");

      auto res = db.insert(config.in_canonical_supercell());
      EXPECT_EQ(res.second, false);

      EXPECT_EQ(config.name(), "SCEL4_2_2_1_1_1_0/0.equiv.0.1");

      EXPECT_EQ(res.first->name(), "SCEL4_2_2_1_1_1_0/0");

      check_made_from_name(config, "SCEL4_2_2_1_1_1_0/0.equiv.0.1");
    }

    {
      // canonical scel, canonical non-primitive occ
      Configuration config(scel);
      config.set_occupation(std::vector<int>({0, 0, 0, 0}));

      EXPECT_EQ(config.name(), "SCEL4_2_2_1_1_1_0/none");

      auto res = db.insert(config.in_canonical_supercell());
      EXPECT_EQ(res.second, true);

      // primitive does not yet exist in database
      EXPECT_EQ(config.name(), "SCEL4_2_2_1_1_1_0/1");
      EXPECT_EQ(res.first->name(), "SCEL4_2_2_1_1_1_0/1");

      // insert primitive
      res = db.insert(config.primitive().in_canonical_supercell());

      // primitive does exist in database
      // - it gets index 2, even though it is the only config from this supercell
      //   in the database, because we don't re-use indices
      EXPECT_EQ(res.first->name(), "SCEL1_1_1_1_0_0_0/2");

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
      config.set_occupation(std::vector<int>({1, 0, 0, 0}));

      // having a different, but equivalent supercell, should not change name ^ see above
      EXPECT_EQ(config.name(), "SCEL4_2_2_1_1_1_0/0");

      auto res = db.insert(config.in_canonical_supercell());
      EXPECT_EQ(res.second, false);
      EXPECT_EQ(config.name(), "SCEL4_2_2_1_1_1_0/0");
      EXPECT_EQ(res.first->name(), "SCEL4_2_2_1_1_1_0/0");

      check_made_from_name(config, "SCEL4_2_2_1_1_1_0/0");
    }

    {
      // non-canonical scel, non-canonical primitive occ
      Configuration config(scel);
      config.set_occupation(std::vector<int>({0, 1, 0, 0}));

      // having a different, but equivalent supercell, should not change name ^ see above
      EXPECT_EQ(config.name(), "SCEL4_2_2_1_1_1_0/0.equiv.0.3");

      auto res = db.insert(config.in_canonical_supercell());
      EXPECT_EQ(res.second, false);

      EXPECT_EQ(config.name(), "SCEL4_2_2_1_1_1_0/0.equiv.0.3");

      EXPECT_EQ(res.first->name(), "SCEL4_2_2_1_1_1_0/0");

      check_made_from_name(config, "SCEL4_2_2_1_1_1_0/0.equiv.0.3");
    }

    {
      // equivalent, but non-canonical scel, canonical non-primitive occ
      Configuration config(scel);
      config.set_occupation(std::vector<int>({0, 0, 0, 0}));

      // having a different, but equivalent supercell, should not change name ^ see above
      EXPECT_EQ(config.name(), "SCEL4_2_2_1_1_1_0/1");

      auto res = db.insert(config.in_canonical_supercell());
      EXPECT_EQ(res.second, false);

      check_made_from_name(config, "SCEL4_2_2_1_1_1_0/super.0.SCEL1_1_1_1_0_0_0/2.equiv.0.0");

    }
  }

  {
    Lattice test_lat = non_canonical_equiv_test_lat(primclex);
    Supercell scel {&primclex, test_lat};

    {
      // non-canonical, non-equivalent, scel, canonical primitive occ
      Configuration config(scel);
      config.set_occupation(std::vector<int>({1, 0, 0, 0, 0, 0, 0, 0}));

      // having a different, but equivalent supercell, should not change name ^ see above
      EXPECT_EQ(config.name(), "SCEL8_4_2_1_1_3_2.4/super.1.SCEL8_4_2_1_1_3_2/none.equiv.0.0");

      auto res = db.insert(config.in_canonical_supercell());
      EXPECT_EQ(res.second, true);
      EXPECT_EQ(config.name(), "SCEL8_4_2_1_1_3_2.4/super.1.SCEL8_4_2_1_1_3_2/0.equiv.0.0");
      EXPECT_EQ(res.first->name(), "SCEL8_4_2_1_1_3_2/0");

      check_made_from_name(config, "SCEL8_4_2_1_1_3_2.4/super.1.SCEL8_4_2_1_1_3_2/0.equiv.0.0");
    }

    {
      // non-canonical, non-equivalent, scel, non-canonical primitive occ
      Configuration config(scel);
      config.set_occupation(std::vector<int>({0, 1, 0, 0, 0, 0, 0, 0}));

      // having a different, but equivalent supercell, should not change name ^ see above
      EXPECT_EQ(config.name(), "SCEL8_4_2_1_1_3_2.4/super.1.SCEL8_4_2_1_1_3_2/0.equiv.0.1");

      auto res = db.insert(config.in_canonical_supercell());
      EXPECT_EQ(res.second, false);

      EXPECT_EQ(config.name(), "SCEL8_4_2_1_1_3_2.4/super.1.SCEL8_4_2_1_1_3_2/0.equiv.0.1");
      EXPECT_EQ(res.first->name(), "SCEL8_4_2_1_1_3_2/0");

      check_made_from_name(config, "SCEL8_4_2_1_1_3_2.4/super.1.SCEL8_4_2_1_1_3_2/0.equiv.0.1");
    }

    {
      // non-canonical, non-equivalent, scel, canonical non-primitive occ
      Configuration config(scel);
      config.set_occupation(std::vector<int>({0, 0, 0, 0, 0, 0, 0, 0}));

      // having a different, but equivalent supercell, should not change name ^ see above
      EXPECT_EQ(config.name(), "SCEL8_4_2_1_1_3_2.4/super.0.SCEL1_1_1_1_0_0_0/2.equiv.0.0");

      auto res = db.insert(config.in_canonical_supercell());
      EXPECT_EQ(res.second, true);

      check_made_from_name(config, "SCEL8_4_2_1_1_3_2.4/super.0.SCEL1_1_1_1_0_0_0/2.equiv.0.0");
    }
  }

}

TEST(ConfigurationTest, Test2) {

  // test Configuration::fill_supercell

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  Supercell scel {&primclex, Lattice(c, a - b, a + b - c)};

  Configuration config(scel);
  EXPECT_EQ(config.size(), 2);

  // include occupation only
  config.set_occupation(std::vector<int>({1, 0}));

  {
    // Identity op
    Supercell scel {&primclex, Lattice(c, a - b, a + b - c)};
    const SymGroup &fg = config.supercell().factor_group();
    Configuration filled = config.fill_supercell(scel, fg[0]);

    Configuration check(scel);
    check.set_occupation(std::vector<int>({1, 0}));

    EXPECT_EQ(filled, check);
  }

  {
    // supercell
    Supercell scel {&primclex, Lattice(c, a - b, 2.*(a + b - c))};
    const SymGroup &fg = config.supercell().factor_group();
    Configuration filled = config.fill_supercell(scel, fg[0]);

    Configuration check(scel);
    check.set_occupation(std::vector<int>({1, 0, 1, 0}));

    EXPECT_EQ(filled, check);
  }

  {
    // 90 deg rotation
    Supercell scel {&primclex, Lattice(c, a - b, a + b - c)};
    const SymGroup &fg = config.supercell().factor_group();
    //std::cout << to_string(fg.info(1), CART) << std::endl;
    Configuration filled = config.fill_supercell(scel, fg[1]);

    Configuration check(scel);
    check.set_occupation(std::vector<int>({1, 0}));

    EXPECT_EQ(filled, check);
  }


  // include occupation & displacements
  Eigen::Vector3d dzero(0., 0., 0.);
  Eigen::Vector3d dx(0.001, 0., 0.);
  Eigen::Vector3d dy(0., 0.001, 0.);
  Eigen::Vector3d dz(0., 0., 0.001);


  //config.set_disp(0, dx);

  {
    // Identity op
    Supercell scel {&primclex, Lattice(c, a - b, a + b - c)};
    const SymGroup &fg = config.supercell().factor_group();
    Configuration filled = config.fill_supercell(scel, fg[0]);

    Configuration check(scel);
    check.set_occupation(std::vector<int>({1, 0}));
    //check.set_disp(0, dx);


    EXPECT_EQ(filled, check);
  }

  {
    // supercell
    Supercell scel {&primclex, Lattice(c, a - b, 2.*(a + b - c))};
    const SymGroup &fg = config.supercell().factor_group();
    Configuration filled = config.fill_supercell(scel, fg[0]);

    Configuration check(scel);
    check.set_occupation(std::vector<int>({1, 0, 1, 0}));
    //check.init_displacement();
    //check.set_disp(0, dx);
    //check.set_disp(2, dx);

    EXPECT_EQ(filled, check);
  }

  {
    // 90 deg rotation
    Supercell scel {&primclex, Lattice(c, a - b, a + b - c)};
    const SymGroup &fg = config.supercell().factor_group();
    Configuration filled = config.fill_supercell(scel, fg[1]);

    Configuration check(scel);
    check.set_occupation(std::vector<int>({1, 0}));
    //check.init_displacement();
    //check.set_disp(0, dy);

    EXPECT_EQ(filled, check);
  }

}

TEST(ConfigurationTest, Test3) {
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
    EXPECT_EQ(scel.is_canonical(), true);

    {
      Configuration config(scel);
      config.set_occupation(std::vector<int>({1, 0, 0, 0}));
      EXPECT_EQ(config.is_canonical(), true);
      EXPECT_EQ(config.is_primitive(), true);
      EXPECT_EQ(config.invariant_subgroup().size(), 48);

      {
        Configuration test(scel);
        test.set_occupation(std::vector<int>({1, 0, 0, 0}));
        EXPECT_EQ(config == test, true);
        EXPECT_EQ(config.is_sym_equivalent(test), true);
        EXPECT_EQ(test.is_sym_equivalent(config), true);
        EXPECT_EQ(test < config, false);
        EXPECT_EQ(config < test, false);
      }

      {
        Configuration test(scel);
        test.set_occupation(std::vector<int>({0, 1, 0, 0}));
        EXPECT_EQ(config == test, false);
        EXPECT_EQ(config.is_sym_equivalent(test), true);
        EXPECT_EQ(test.is_sym_equivalent(config), true);
        EXPECT_EQ(test < config, true);
        EXPECT_EQ(config < test, false);

        auto to_canonical = test.to_canonical();
        EXPECT_EQ(to_canonical.factor_group_index(), 0);
        EXPECT_EQ(to_canonical.translation_index(), 1);
        EXPECT_EQ(copy_apply(to_canonical, test) == config, true);
      }
    }
  }
}
