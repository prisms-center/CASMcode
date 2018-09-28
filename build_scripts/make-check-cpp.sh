# configure && make && make check, with helpful messages
#
# Some tests are only run if certains conditions exist:
# - Tests against test projects are run if CASM_TEST_PROJECTS_DIR exists
# - Some vasp tests are run if a 'vasp' executable is found

set -e
detect_os

check_var "CASM_BUILD_DIR" "CASMcode repository location"
check_var "CASM_TESTS" "Particular test categories to run (default="", runs all tests)" ""
check_var "CASM_TEST_FLAGS" "Customize the options given to the test programs" "--log_level=test_suite --catch_system_errors=no"

. $CASM_BUILD_DIR/build_scripts/make-cpp.sh

# Run tests and print output
cd $CASM_BUILD_DIR
echo "CASM_TESTS: '$CASM_TESTS'"
echo "CASM_TEST_FLAGS: '$CASM_TEST_FLAGS'"
echo "CASM_BOOST_PREFIX: '$CASM_BOOST_PREFIX'"
echo 'make check -j $CASM_NCPU CASM_BOOST_PREFIX=$CASM_BOOST_PREFIX TESTS="$CASM_TESTS" TEST_FLAGS="$CASM_TEST_FLAGS"'

rm -rf test-suite.log
make check -j $CASM_NCPU CASM_BOOST_PREFIX="$CASM_BOOST_PREFIX" ${CASM_TESTS:+TESTS="$CASM_TESTS"} TEST_FLAGS="$CASM_TEST_FLAGS" \
  || do_if_failed "eval [ -f test-suite.log ] && cat test-suite.log" \
  || { echo "'make check -j $CASM_NCPU CASM_BOOST_PREFIX=$CASM_BOOST_PREFIX' failed"; exit 1; }
