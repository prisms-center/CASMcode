# for running "make check" or "pytest" for development in a conda environment
# - the conda environment must be activated

### initialization - shouldn't need to touch
set -e
export CASM_BUILD_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
. $CASM_BUILD_DIR/build_scripts/install-functions.sh
detect_os
check_var "CONDA_PREFIX" "Must have the conda environment activated"
export CASM_PREFIX=$CONDA_PREFIX

### end initialization ###

### variables - Control how CASM is built and what is tested ###

# For C++ and Python
check_var "CASM_TEST_PROJECTS_DIR" "Location of CASM_test_projects" ""
check_var "CASM_SKIP_MAKE" "Set to non-zero length to skip initialize make" ""

# For C++ tests
check_var "CASM_SKIP_CPP_TESTS" "Set to non-zero length to skip ++ tests" ""
check_var "CASM_CXXFLAGS" "Compiler flags" ""
check_var "CASM_NCPU" "Compiler -j option" 2
check_var "CASM_TESTS" "Particular test categories to run (default="", runs all tests)" ""
check_var "CASM_TEST_FLAGS" "Customize the options given to the test programs" "--log_level=test_suite --catch_system_errors=no"

# For Python tests
check_var "CASM_SKIP_PYTHON_TESTS" "Set to non-zero length to skip Python tests" ""
check_var "CASM_PYTEST_ARGS" "Arguments to pass to pytest" "test_casm"

### end variables ###

# set OS-dependent variable defaults
#   only CASM_CONFIGFLAGS can't be overridden from this script
. $CASM_BUILD_DIR/build_scripts/travis-variables-$CASM_OS_NAME.sh

# Set CASM_SKIP_MAKE to skip make
if ! [ -n "$CASM_SKIP_MAKE" ]; then
  bash $CASM_BUILD_DIR/build_scripts/make-cpp.sh
else
  echo "skipping make"
fi

# Set CASM_SKIP_CPP to skip C++ tests
if ! [ -n "$CASM_SKIP_CPP_TESTS" ]; then
  bash $CASM_BUILD_DIR/build_scripts/make-check-cpp.sh
else
  echo "skipping c++ tests"
fi

# Set CASM_SKIP_PYTHON to skip Python tests (Cpp must be built to run Python tests)
if ! [ -n "$CASM_SKIP_PYTHON_TESTS" ]; then
  bash $CASM_BUILD_DIR/build_scripts/check-python.sh "${PYTEST_ARGS[@]}"
else
  echo "skipping Python tests"
fi
