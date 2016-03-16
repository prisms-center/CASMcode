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
  
  std::vector<int> sublat_indices = {0};
  PrimNeighborList::Matrix3Type W;
  W.row(0) << 2, 1, 1;
  W.row(1) << 1, 2, 1;
  W.row(2) << 1, 1, 2;
  
  PrimNeighborList nlist(W, sublat_indices.begin(), sublat_indices.end());
  
  Clexulator clexulator("test_Clexulator",
                        "tests/unit/clex",
                        nlist,
                        compile_opt,
                        so_opt);

  BOOST_CHECK_EQUAL(clexulator.corr_size(), 75);

}

BOOST_AUTO_TEST_SUITE_END()
