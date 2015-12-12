#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/Clexulator.hh"

/// Dependencies

/// What is being used to test it:
#include <boost/filesystem.hpp>

using namespace CASM;

BOOST_AUTO_TEST_SUITE(ClexulatorTest)

BOOST_AUTO_TEST_CASE(MakeClexulatorTest) {
  namespace fs = boost::filesystem;
  
  std::string compile_opt = RuntimeLibrary::cxx() + " " + RuntimeLibrary::default_cxxflags() + " -Iinclude";
  std::string so_opt = RuntimeLibrary::default_so_options() + " -lboost_system";

  if(std::getenv("CASMBOOST_PATH") != nullptr) {
    fs::path boost_path(std::getenv("CASMBOOST_PATH"));
    compile_opt += " -I" + (boost_path / "include").string();
    so_opt += " -L" + (boost_path / "lib").string();
  }
  
  Clexulator clexulator("test_Clexulator",
                        "tests/unit/clex",
                        compile_opt,
                        so_opt);

  BOOST_CHECK_EQUAL(clexulator.corr_size(), 75);

}

BOOST_AUTO_TEST_SUITE_END()
