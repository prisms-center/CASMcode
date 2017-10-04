CASM tests
==========

`CASM` uses the [boost unit test framework](http://www.boost.org/doc/libs/1_42_0/libs/test/doc/html/index.html) for testing.


Testing dependencies
--------------------

No differences from installation.


Testing configuration
---------------------

No differences from installation.

Writing Tests
-------------

See unit tests examples in ``tests/unit/<category>/test_<name>.cpp``:

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

BOOST_AUTO_TEST_CASE(Test1) {

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
    make check
```

Run particular test categories
---------------------------------------
- To run particular test categories, set the ``TESTS`` variable

From ``CASMcode`` directory:

```
	TEST1=tests/unit/App/run_test_App
	TEST2=tests/unit/clex/run_test_clex
	...
	make check TESTS="$TEST1 $TEST2 ..."
```

Run particular test cases
-------------------------

- To customize the options given to the test programs, for instance to run only particular test suites or test cases, set ``TEST_FLAGS``.
- For more options see [docs](http://www.boost.org/doc/libs/1_42_0/libs/test/doc/html/utf/user-guide/runtime-config/run-by-name.html)
 
From ``CASMcode`` directory:

```
    TEST1=tests/unit/clex/run_test_clex
    TEST_FLAGS="--run_test=ScelEnumEquivalentsTest/Test2 "
    TEST_FLAGS+="--log_level=test_suite "
    TEST_FLAGS+="--catch_system_errors=no "
    make check TESTS="$TEST1" TEST_FLAGS="$TEST_FLAGS"
```

The default ``TEST_FLAGS`` are:

```
    TEST_FLAGS="--log_level=test_suite --catch_system_errors=no"
```

