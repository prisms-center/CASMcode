#include "gtest/gtest.h"
#include "autotools.hh"

/// What is being tested:
#include "casm/system/RuntimeLibrary.hh"

/// What is being used to test it:
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

using namespace CASM;

TEST(RuntimeLibraryTest, FunctionTest) {

  EXPECT_EQ(true, true);
  std::string cc_filename_base = autotools::abs_srcdir() + "/tests/unit/system/runtime_lib";
  fs::path cc_filename {cc_filename_base + ".cc"};
  EXPECT_EQ(true, true);

  fs::ofstream file(cc_filename);
  file << "#include <iostream>\n"
       "extern \"C\" int forty_two() {\n"
       "   return 42;\n"
       "}\n"
       "\n"
       "extern \"C\" int add(int a, int b) {\n"
       "   return a + b;\n"
       "}\n";
  file.close();
  EXPECT_EQ(true, true);

  std::string compile_opt = RuntimeLibrary::default_cxx().first + " " +
                            RuntimeLibrary::default_cxxflags().first;
  std::string so_opt = RuntimeLibrary::default_cxx().first + " " +
                       RuntimeLibrary::default_soflags().first + " " +
                       link_path(RuntimeLibrary::default_boost_libdir().first.string());
  EXPECT_EQ(true, true);

  RuntimeLibrary lib(
    cc_filename_base,
    compile_opt,
    so_opt);

  EXPECT_EQ(true, true);

  // get the 'int forty_two()' function
  std::function<int()> forty_two = lib.get_function<int()>("forty_two");

  // use it to do something
  EXPECT_EQ(42, forty_two());

  // get the 'int add(int, int)' function
  std::function<int(int, int)> add = lib.get_function<int(int, int)>("add");

  // use it to do something
  EXPECT_EQ(5, add(2, 3));

  // delete the library
  lib.rm();

  EXPECT_EQ(true, true);

}
