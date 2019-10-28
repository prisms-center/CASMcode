#include "gtest/gtest.h"

/// What is being tested:
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/ScelEnum_impl.hh"

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

  test::ZrOProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  primclex.settings().set_crystallography_tol(1e-5);

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  std::vector<std::string> m_names;

  // -- Test ScelEnumByProps --------------------
  {
    ScelEnumProps enum_props(1, 10);
    ScelEnumByProps e(primclex, enum_props);

    EXPECT_EQ(e.name(), "ScelEnumByProps");

    auto it = e.begin();
    EXPECT_TRUE(true);
    EXPECT_EQ(it.name(), "ScelEnumByProps");

    auto end = e.end();
    EXPECT_TRUE(true);

    Index count = 0;
    for(; it != end; ++it, ++count) {
      m_names.push_back(it->name());

      Lattice canon_check = canonical_equivalent_lattice(
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

  // -- use results to Test ScelEnumByName --------------------
  {
    ScelEnumByName e(primclex, m_names.begin(), m_names.end());
    EXPECT_EQ(e.name(), "ScelEnumByName");

    auto it = e.begin();
    EXPECT_TRUE(true);
    EXPECT_EQ(it.name(), "ScelEnumByName");

    auto end = e.end();
    EXPECT_TRUE(true);

    Index count = 0;
    for(; it != end; ++it, ++count) {
      //std::cout << it->name() << std::endl;
    }
    EXPECT_EQ(count, 114);
    EXPECT_TRUE(it == end);
  }

}

TEST(ScelEnumTest, Test2) {

  // create a project
  test::FCCTernaryProj proj;
  proj.check_init();

  // in case you want to see what's happening
  OStringStreamLog ss_log;
  OStringStreamLog ss_debug_log;
  OStringStreamLog ss_err_log;

  // construct PrimClex
  PrimClex primclex(proj.dir, Logging(ss_log, ss_debug_log, ss_err_log));

  auto exec = [&](const std::string & args) {
    CommandArgs cmdargs(args, &primclex, proj.dir, ss_log, ss_err_log);
    int code = casm_api(cmdargs);
    //std::cout << ss_log.ss().str() << std::endl;
    //std::cout << ss_err_log.ss().str() << std::endl;
    return code;
  };

  EXPECT_EQ(exec("casm enum -h"), 0);
  EXPECT_EQ(exec("casm enum --method ScelEnum --max 4"), 0);
  EXPECT_EQ(exec("casm enum --method ConfigEnumAllOccupations --all"), 0);
  EXPECT_EQ(exec("casm enum --method ScelEnum --max 8"), 0);
  EXPECT_EQ(exec("casm enum --method ConfigEnumAllOccupations --max 6 -i '{\"existing_only\":true}'"), 0);

}
