#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/system/RuntimeLibrary.hh"

/// Dependencies

/// What is being used to test it:
#include <boost/filesystem.hpp>

using namespace CASM;

BOOST_AUTO_TEST_SUITE(RuntimeLibraryTest)

BOOST_AUTO_TEST_CASE(FunctionTest) {
  
  std::string cc_filename_base = "runtime_lib";

  std::string cc_file;
  cc_file = std::string("#include <iostream>\n") +
            "extern \"C\" int forty_two() {\n" +
            "   return 42;\n" +
            "}\n" +
            "\n" +
            "extern \"C\" int add(int a, int b) {\n" +
            "   return a + b;\n" +
            "}\n";

  RuntimeLibrary lib;

  // write the library file and compile
  lib.compile("tests/unit/system/runtime_lib", cc_file.c_str());

  // load the library
  lib.load("tests/unit/system/runtime_lib");

  // get the 'int forty_two()' function
  std::function<int()> forty_two = lib.get_function<int()>("forty_two");

  // use it to do something
  BOOST_CHECK_EQUAL(42, forty_two());
  
  // get the 'int add(int, int)' function
  std::function<int(int, int)> add = lib.get_function<int(int, int)>("add");

  // use it to do something
  BOOST_CHECK_EQUAL(5, add(2,3));
  
  // close the library
  lib.close();

  // delete the library
  lib.rm();
  
}

BOOST_AUTO_TEST_SUITE_END()
