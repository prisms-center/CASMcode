#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/system/RuntimeLibrary.hh"

/// Dependencies
#include "casm/casm_io/Log.hh"

/// What is being used to test it:
#include <boost/filesystem.hpp>

using namespace CASM;

BOOST_AUTO_TEST_SUITE(RuntimeLibraryTest)

BOOST_AUTO_TEST_CASE(FunctionTest) {

  BOOST_CHECK_EQUAL(true, true);
  std::string cc_filename_base = "tests/unit/system/runtime_lib";
  fs::path cc_filename {cc_filename_base + ".cc"};
  BOOST_CHECK_EQUAL(true, true);

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
  BOOST_CHECK_EQUAL(true, true);

  std::string compile_opt = RuntimeLibrary::default_cxx().first + " " +
                            RuntimeLibrary::default_cxxflags().first;
  std::string so_opt = RuntimeLibrary::default_cxx().first + " " +
                       RuntimeLibrary::default_soflags().first + " " +
                       link_path(RuntimeLibrary::default_boost_libdir().first.string());
  BOOST_CHECK_EQUAL(true, true);

  RuntimeLibrary lib(
    cc_filename_base,
    compile_opt,
    so_opt,
    "Compiling RuntimeLibrary test code",
    Logging());

  BOOST_CHECK_EQUAL(true, true);

  // get the 'int forty_two()' function
  std::function<int()> forty_two = lib.get_function<int()>("forty_two");

  // use it to do something
  BOOST_CHECK_EQUAL(42, forty_two());

  // get the 'int add(int, int)' function
  std::function<int(int, int)> add = lib.get_function<int(int, int)>("add");

  // use it to do something
  BOOST_CHECK_EQUAL(5, add(2, 3));

  // delete the library
  lib.rm();

  BOOST_CHECK_EQUAL(true, true);

}

BOOST_AUTO_TEST_SUITE_END()
