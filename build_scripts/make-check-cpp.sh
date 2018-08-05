# configure && make && make check, with helpful messages
#
# Some tests are only run if certains conditions exist:
# - Tests against test projects are run if CASM_TEST_PROJECTS_DIR exists
# - Some vasp tests are run if a 'vasp' executable is found

set -e
detect_os

check_var "CASM_BUILD_DIR" "CASMcode repository location"

. $CASM_BUILD_DIR/build_scripts/make-cpp.sh

# Run tests and print output
cd $CASM_BUILD_DIR
make check -j $CASM_NCPU CASM_BOOST_PREFIX="$CASM_BOOST_PREFIX" \
  || bash ./build_scripts/check-rpath.sh \
  || do_if_failed "cat $CASM_BUILD_DIR/test-suite.log" \
  || { echo "'make check' failed"; exit 1; }
