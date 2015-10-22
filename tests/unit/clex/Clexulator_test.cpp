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
  
  std::string boost_path = "";
  if(std::getenv("CASMBOOST_PATH") != nullptr) {
    boost_path = (fs::path(std::getenv("CASMBOOST_PATH")) / "lib").string();
  }
  
  std::string compile_opt = RuntimeLibrary::default_compile_options() + " -Iinclude";
  std::string so_opt = RuntimeLibrary::default_so_options() + " -lboost_system -L" + boost_path;
  
  Clexulator clexulator("test_Clexulator",
                        "tests/unit/clex",
                        compile_opt,
                        so_opt);

  BOOST_CHECK_EQUAL(clexulator.corr_size(), 75);

}

BOOST_AUTO_TEST_SUITE_END()
