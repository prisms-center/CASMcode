# use CASM conda development environment to build and run all tests for casm-cpp and casm-python
#   if CASM_TEST_PROJECTS_DIR is non-zero length, requires
#     CASM_TEST_PROJECTS_ID and MC_API_KEY

set -e
export CASM_BUILD_DIR=${CASM_BUILD_DIR:-${TRAVIS_BUILD_DIR:-$(pwd)}}
. $CASM_BUILD_DIR/build_scripts/install-functions.sh
. $CASM_BUILD_DIR/build_scripts/build_versions.sh
detect_os

check_var "CASM_CONDA_DIR" "Location to install conda and conda environments" "$HOME/.local/conda"
check_var "CASM_VERSION" "CASM version (used for naming conda env)" "0.3"
check_var "CASM_PYTHON_VERSION" "Python version"
check_var "CASM_ENV_NAME" "CASM conda environment name" "casm_${CASM_VERSION}_py${CASM_PYTHON_VERSION}"
check_var "CASM_TEST_PROJECTS_DIR" "Location to download CASM_test_projects" ""

# activate conda
. $CASM_CONDA_DIR/etc/profile.d/conda.sh
conda activate $CASM_ENV_NAME

export CASM_PREFIX=$CONDA_PREFIX

# set OS-dependent variables
. $CASM_BUILD_DIR/build_scripts/travis-variables-$CASM_OS_NAME.sh

# make-check-cpp
bash $CASM_BUILD_DIR/build_scripts/make-check-cpp.sh \
  || { echo "'make-check-cpp.sh' failed"; exit 1; }

# check-python
bash $CASM_BUILD_DIR/build_scripts/check-python.sh \
  || { echo "'check-python.sh' failed"; exit 1; }
