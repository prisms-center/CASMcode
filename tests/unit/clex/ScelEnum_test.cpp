#include "gtest/gtest.h"

/// What is being tested:
#include "casm/clex/ScelEnum.hh"

/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "ZrOProj.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/app/casm_functions.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/Niggli.hh"

using namespace CASM;

TEST(ScelEnumTest, Test1) {

  ScopedNullLogging logging;
  test::ZrOProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir);
  primclex.settings().set_crystallography_tol(1e-5);

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  // -- Test ScelEnumByProps --------------------
  {
    xtal::ScelEnumProps enum_props(1, 10);
    ScelEnumByProps e(primclex.shared_prim(), enum_props);

    EXPECT_EQ(e.name(), "ScelEnumByProps");

    auto it = e.begin();
    EXPECT_TRUE(true);
    EXPECT_EQ(it.name(), "ScelEnumByProps");

    auto end = e.end();
    EXPECT_TRUE(true);

    Index count = 0;
    for(; it != end; ++it, ++count) {

      Lattice canon_check = xtal::canonical::equivalent(
                              it->lattice(),
                              primclex.prim().point_group(),
                              primclex.crystallography_tol());

      bool check = almost_equal(
                     it->lattice().lat_column_mat(),
                     canon_check.lat_column_mat(),
                     primclex.crystallography_tol());

      if(!check) {
        std::cout << "superlat: \n" << it->lattice().lat_column_mat() << std::endl;
        std::cout << "canon_check: \n" << canon_check.lat_column_mat() << std::endl;
      }

      EXPECT_EQ(check, true);

      EXPECT_EQ(it->is_canonical(), true);

    }
    EXPECT_EQ(count, 114);
    EXPECT_TRUE(it == end);
  }

}

TEST(ScelEnumTest, Test2) {

  // in case you want to see what's happening
  ScopedStringStreamLogging logging;

  // create a project
  test::FCCTernaryProj proj;
  proj.check_init();

  // construct PrimClex
  PrimClex primclex(proj.dir);

  auto exec = [&](const std::string & args) {
    // std::cout << "\n---------------\n" << std::endl;
    // std::cout << "args: " << args << std::endl;
    CommandArgs cmdargs(args, &primclex, proj.dir);
    int code = casm_api(cmdargs);
    // std::cout << "\n---------------\n" << std::endl;
    // std::cout << "log: " << std::endl;
    // std::cout << logging.ss().str() << std::endl;
    // std::cout << "\n---------------\n" << std::endl;
    // std::cout << "err_log: " << std::endl;
    // std::cout << logging.err_ss().str() << std::endl;
    return code;
  };

  EXPECT_EQ(exec("casm enum -h"), 0);
  EXPECT_EQ(exec("casm enum --method ScelEnum --max 4"), 0);
  EXPECT_EQ(exec("casm enum --method ConfigEnumAllOccupations --all"), 0);
  EXPECT_EQ(exec("casm enum --method ScelEnum --max 8"), 0);
  EXPECT_EQ(exec("casm enum --method ConfigEnumAllOccupations --max 6"), 0);

}
