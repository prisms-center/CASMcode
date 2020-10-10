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
  ScopedStringStreamLogging logging;

  // construct PrimClex
  PrimClex primclex(proj.dir);

  auto exec = [&](const std::string & args) {
    CommandArgs cmdargs(args, &primclex, proj.dir);
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

  //std::cout << logging.ss().str() << std::endl;
  //std::cout << logging.err_ss().str() << std::endl;

}
