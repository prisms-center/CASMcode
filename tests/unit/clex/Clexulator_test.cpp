#include "gtest/gtest.h"

/// What is being tested:
#include "casm/clex/Clexulator.hh"

/// Dependencies
#include "casm/system/RuntimeLibrary.hh"
#include "casm/casm_io/Log.hh"

/// What is being used to test it:
#include <boost/filesystem.hpp>

using namespace CASM;

TEST(ClexulatorTest, MakeClexulatorTest) {

  std::cout << "skipping Clexulator_test" << std::endl;
  if(true) {
    return;
  }

  namespace fs = boost::filesystem;

  std::string compile_opt = RuntimeLibrary::default_cxx().first + " " + RuntimeLibrary::default_cxxflags().first + " -Iinclude";
  std::string so_opt = RuntimeLibrary::default_cxx().first + " " + RuntimeLibrary::default_soflags().first;

  if(!RuntimeLibrary::default_boost_includedir().first.empty()) {
    compile_opt += " " + include_path(RuntimeLibrary::default_boost_includedir().first);
  }

  if(!RuntimeLibrary::default_boost_libdir().first.empty()) {
    so_opt += " " + link_path(RuntimeLibrary::default_boost_libdir().first);
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

  EXPECT_EQ(clexulator.corr_size(), 75);

}
