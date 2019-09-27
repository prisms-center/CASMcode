#include "gtest/gtest.h"

/// What is being tested:
#include "casm/app/casm_functions.hh"

/// What is being used to test it:

//#include "casm/app/ProjectBuilder.hh"
#include "Common.hh"
#include "FCCTernaryProj.hh"

using namespace CASM;

TEST(settingsTest, Basics) {

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
    return casm_api(cmdargs);
  };

  EXPECT_EQ(exec("casm settings -l"), 0);
  EXPECT_EQ(exec("casm settings --new-calctype test1"), 0);
  EXPECT_EQ(exec("casm settings --new-ref test1"), 0);
  EXPECT_EQ(exec("casm settings --new-bset test1"), 0);
  EXPECT_EQ(exec("casm settings --new-eci test1"), 0);
  EXPECT_EQ(exec("casm settings --set-calctype default"), 0);
  EXPECT_EQ(exec("casm settings --set-calctype does_not_exist"), 1);
  EXPECT_EQ(exec("casm settings --set-calctype test1"), 0);
  EXPECT_EQ(exec("casm settings --set-bset test1"), 0);
  EXPECT_EQ(exec("casm settings --set-bset default"), 0);

  //std::cout << ss_log.ss().str() << std::endl;
  //std::cout << ss_err_log.ss().str() << std::endl;

}
