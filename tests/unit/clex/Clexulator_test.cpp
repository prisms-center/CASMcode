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

  std::string compile_opt = RuntimeLibrary::default_cxx().first + " " + RuntimeLibrary::default_cxxflags().first + " -Iinclude";
  std::string so_opt = RuntimeLibrary::default_cxx().first + " " + RuntimeLibrary::default_soflags().first;

  if(!RuntimeLibrary::default_boost_prefix().first.empty()) {
    fs::path boost_prefix = RuntimeLibrary::default_boost_prefix().first;
    compile_opt += " " + include_path(boost_prefix);
    so_opt += " " + link_path(boost_prefix);
  }

  std::vector<int> sublat_indices = {0};
  PrimNeighborList::Matrix3Type W;
  W.row(0) << 2, 1, 1;
  W.row(1) << 1, 2, 1;
  W.row(2) << 1, 1, 2;

  PrimNeighborList nlist(W, sublat_indices.begin(), sublat_indices.end());

  Log dumblog = null_log();

  Clexulator clexulator("test_Clexulator",
                        "tests/unit/clex",
                        nlist,
                        dumblog,
                        compile_opt,
                        so_opt);

  BOOST_CHECK_EQUAL(clexulator.corr_size(), 75);

}

BOOST_AUTO_TEST_SUITE_END()
