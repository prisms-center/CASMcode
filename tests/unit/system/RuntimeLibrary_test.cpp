#include "casm/system/RuntimeLibrary.hh"

#include <filesystem>
#include <fstream>

#include "Common.hh"
#include "gtest/gtest.h"

using namespace CASM;

TEST(RuntimeLibraryTest, FunctionTest) {
  EXPECT_EQ(true, true);
  test::TmpDir tmpdir;
  std::string cc_filename_base = (tmpdir.path() / "runtime_lib").string();
  fs::path cc_filename{cc_filename_base + ".cc"};
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
  EXPECT_TRUE(fs::exists(cc_filename)) << "does not exist: " << cc_filename;

  // std::string compile_opt = RuntimeLibrary::default_cxx().first + " " +
  //                           RuntimeLibrary::default_cxxflags().first;
  std::string compile_opt =
      RuntimeLibrary::default_cxx().first + " -O3 -Wall -fPIC --std=c++17 ";
  std::cout << "compile_opt: " << compile_opt << std::endl;

  // std::string so_opt =
  //     RuntimeLibrary::default_cxx().first + " " +
  //     RuntimeLibrary::default_soflags().first + " " +
  //     link_path(RuntimeLibrary::default_boost_libdir().first.string());
  std::string so_opt = RuntimeLibrary::default_cxx().first + " -shared ";
  std::cout << "so_opt: " << so_opt << std::endl;
  EXPECT_EQ(true, true);

  try {
    RuntimeLibrary lib(cc_filename_base, compile_opt, so_opt);

    EXPECT_EQ(true, true);

    // get the 'int forty_two()' function
    std::function<int()> forty_two = lib.get_function<int()>("forty_two");

    // use it to do something
    EXPECT_EQ(42, forty_two());

    // get the 'int add(int, int)' function
    std::function<int(int, int)> add = lib.get_function<int(int, int)>("add");

    // use it to do something
    EXPECT_EQ(5, add(2, 3));

    // TODO: This causes googletest to hit a segmentation fault when exiting the
    // test case
    // delete the library
    lib.rm();

    EXPECT_FALSE(fs::exists(cc_filename)) << "does exist: " << cc_filename;

  } catch (runtime_lib_compile_error &e) {
    e.print(std::cout);
    throw;
  } catch (runtime_lib_shared_error &e) {
    e.print(std::cout);
    throw;
  } catch (std::exception &e) {
    std::cout << "error: " << e.what() << std::endl;
    throw;
  }
}
