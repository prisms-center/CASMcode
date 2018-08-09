CASM C++ tests
==========

`CASM` uses the [boost unit test framework](http://www.boost.org/doc/libs/1_42_0/libs/test/doc/html/index.html) for C++ testing.


Writing Tests
-------------

See unit tests examples in ``tests/unit/<category>/<name>_test.cpp``:

- Boost unit test framework [documentation](http://www.boost.org/doc/libs/1_42_0/libs/test/doc/html/index.html)
- Tests are organized by category into directories matching the CASM source code. Within a category there is typically one test suite per file with a name matching the name used for header and source files of the code being tested.
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


Run particular test categories
---------------------------------------
- To run particular C++ test categories, set the ``CASM_TESTS`` environment variable

From ``CASMcode`` directory:

```
$ TEST1=tests/unit/App/run_test_App
$ TEST2=tests/unit/clex/run_test_clex
  ...
$ export CASM_TESTS="$TEST1 $TEST2 ..."
$ bash build_test.sh
```

Run particular test cases
-------------------------

- To customize the options given to the test programs, for instance to run only particular test suites or test cases, set the ``CASM_TEST_FLAGS`` environment variable.
- For more options see [docs](http://www.boost.org/doc/libs/1_42_0/libs/test/doc/html/utf/user-guide/runtime-config/run-by-name.html)
 
From ``CASMcode`` directory:

```
$ export CASM_TESTS="tests/unit/clex/run_test_clex"
$ CASM_TEST_FLAGS="--run_test=ScelEnumEquivalentsTest/Test2 "
$ CASM_TEST_FLAGS+="--log_level=test_suite "
$ CASM_TEST_FLAGS+="--catch_system_errors=no "
$ export CASM_TEST_FLAGS
$ 
$ bash build_test.sh
```

The default ``CASM_TEST_FLAGS`` are:

```
$ CASM_TEST_FLAGS="--log_level=test_suite --catch_system_errors=no"
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