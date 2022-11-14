#include "gtest/gtest.h"

/// What is being tested:
#include "casm/app/LogRuntimeLibrary.hh"
#include "casm/system/RuntimeLibrary.hh"

/// Dependencies
#include "Common.hh"
#include "casm/casm_io/Log.hh"

using namespace CASM;

TEST(LogRuntimeLibraryTest, FunctionTest) {
  ScopedStringStreamLogging logging;

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
  EXPECT_EQ(true, true);

  std::string compile_opt =
      RuntimeLibrary::default_cxx().first + " -O3 -Wall -fPIC --std=c++17 ";
  std::cout << "compile_opt: " << compile_opt << std::endl;

  std::string so_opt = RuntimeLibrary::default_cxx().first + " -shared ";
  std::cout << "so_opt: " << so_opt << std::endl;
  EXPECT_EQ(true, true);

  std::shared_ptr<RuntimeLibrary> lib =
      log_make_shared_runtime_lib(cc_filename_base, compile_opt, so_opt,
                                  "Compiling RuntimeLibrary test code");

  EXPECT_EQ(true, true);

  // get the 'int forty_two()' function
  std::function<int()> forty_two = lib->get_function<int()>("forty_two");

  // use it to do something
  EXPECT_EQ(42, forty_two());

  // get the 'int add(int, int)' function
  std::function<int(int, int)> add = lib->get_function<int(int, int)>("add");

  // use it to do something
  EXPECT_EQ(5, add(2, 3));

  // delete the library
  lib->rm();

  EXPECT_EQ(true, true);
}
