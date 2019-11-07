CASM C++ tests
==========

`CASM` uses the [googletest](https://github.com/google/googletest/blob/master/googletest/docs/primer.md) testing framework for C++ testing.

Selecting tests: https://github.com/google/googletest/blob/master/googletest/docs/advanced.md#selecting-tests


Writing Tests
-------------

See unit tests examples in ``tests/unit/<group>/<name>_test.cpp``:

- Boost unit test framework [documentation](http://www.boost.org/doc/libs/1_42_0/libs/test/doc/html/index.html)
- Tests are grouped into directories matching the CASM source code. Within a group there is typically one test suite per file with a name matching the name used for header and source files of the code being tested.
- Test projects, either ZrO (``test::ZrOProj``) or a generic FCC ternary (``test::FCCTernary``),  and their associated ``PrimClex`` can be constructed to use for testing.
- An example test using the ZrO test project:

```
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Common.hh"
#include "ZrOProj.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(MyClassTest)

BOOST_AUTO_TEST_CASE(Test1_ZrO) {

    test::ZrOProj proj;
    proj.check_init();
    proj.check_composition(); # if composition axes required
    PrimClex primclex(proj.dir, null_log());
    ...
}

BOOST_AUTO_TEST_SUITE_END()
```


Run all tests
-------------

From ``CASMcode`` directory:

```
$ bash build_test.sh
```

Options:

- Set ``CASM_SKIP_CPP_TESTS`` to non-zero length to skip C++ tests.
- Do ``unset CASM_SKIP_CPP_TESTS`` to re-enable C++ tests.
- Set ``CASM_SKIP_PYTHON_TESTS`` to non-zero length to skip Python tests.
- Do ``unset CASM_SKIP_PYTHON_TESTS`` to re-enable Python tests.
- Use the ``CASM_PYTEST_ARGS`` environment variable to pass arguments to pytest.
- Use the ``CASM_PYTEST_OPTIONS`` environment variable to set pytest options. Default is ``"-r ap -s"``
- Add a file named `'skip'` to a test directory to skip all test classes in that directory.
- Add a file named `'skip_TestCase'` to skip the `TestCase` class in that directory.


Run particular test programs
----------------------------
- CASM tests are grouped into separate test programs for different CASM components
- To run a particular C++ test program, set the ``CASM_TESTS`` environment variable

From ``CASMcode`` directory:

```
$ TEST1=casm_unit_crystallography
$ TEST2=casm_unit_symmetry
  ...
$ export CASM_TESTS="$TEST1 $TEST2 ..."
$ bash build_test.sh
```

Run particular test cases
-------------------------

- To customize the options given to the test programs, for instance to run only particular test suites or test cases, set the ``CASM_TEST_FLAGS`` environment variable.
- For more options see [selecting tests](https://github.com/google/googletest/blob/master/googletest/docs/primer.md)


From ``CASMcode`` directory:

```
$ export CASM_TESTS="casm_test_crystallography "
$ export CASM_TEST_FLAGS="--gtest_filter=CoordinateTest.* --gtest_color=yes "
$
$ bash build_test.sh
```

Clean test output
-----------------
From ``CASMcode`` directory:

```
$ bash checkclean.sh
```

CI Testing
----------
Latest release: [![Build Status](https://travis-ci.org/prisms-center/CASMcode.svg?branch=master)](https://travis-ci.org/prisms-center/CASMcode)

Development: [![Build Status](https://travis-ci.org/prisms-center/CASMcode.svg?branch=0.2.X)](https://travis-ci.org/prisms-center/CASMcode)
