#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/app/casm_functions.hh"

/// What is being used to test it:

//#include "casm/app/ProjectBuilder.hh"
#include "Common.hh"
#include "FCCTernaryProj.hh"

using namespace CASM;


BOOST_AUTO_TEST_SUITE(settingsTest)

BOOST_AUTO_TEST_CASE(Basics) {

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

  BOOST_CHECK_EQUAL(exec("casm settings -l"), 0);
  BOOST_CHECK_EQUAL(exec("casm settings --new-calctype test1"), 0);
  BOOST_CHECK_EQUAL(exec("casm settings --new-ref test1"), 0);
  BOOST_CHECK_EQUAL(exec("casm settings --new-bset test1"), 0);
  BOOST_CHECK_EQUAL(exec("casm settings --new-eci test1"), 0);
  BOOST_CHECK_EQUAL(exec("casm settings --set-calctype default"), 0);
  BOOST_CHECK_EQUAL(exec("casm settings --set-calctype does_not_exist"), 1);
  BOOST_CHECK_EQUAL(exec("casm settings --set-calctype test1"), 0);
  BOOST_CHECK_EQUAL(exec("casm settings --set-bset test1"), 0);
  BOOST_CHECK_EQUAL(exec("casm settings --set-bset default"), 0);

  //std::cout << ss_log.ss().str() << std::endl;
  //std::cout << ss_err_log.ss().str() << std::endl;

}

BOOST_AUTO_TEST_SUITE_END()