# configure && make && make check, with helpful messages
#
# Some tests are only run if certains conditions exist:
# - Some vasp tests are run if a 'vasp' executable is found

set -e
detect_os

check_var "CASM_BUILD_DIR" "CASMcode repository location"
check_var "CASM_TESTS" "Particular test categories to run (default="", runs all tests)" ""

. $CASM_BUILD_DIR/build_scripts/make-cpp.sh

# Run tests and print output
cd $CASM_BUILD_DIR
echo "CASM_TESTS: '$CASM_TESTS'"
echo 'make check -j $CASM_NCPU ${CASM_TESTS:+TESTS="$CASM_TESTS"}'

rm -rf test-suite.log
make check -j $CASM_NCPU ${CASM_TESTS:+TESTS="$CASM_TESTS"} \
  || do_if_failed "eval [ -f test-suite.log ] && cat test-suite.log" \
  || { echo "'make check -j $CASM_NCPU ${CASM_TESTS:+TESTS="$CASM_TESTS"}' failed"; exit 1; }
