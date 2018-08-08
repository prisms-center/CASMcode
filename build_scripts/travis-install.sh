# setup environment and run conda-devel.sh
#   should work on local machine (run from CASMcode repo or with CASM_BUILD_DIR set) or travis for osx and linux
#   if CASM_TEST_PROJECTS_DIR is non-zero length, requires
#     CASM_TEST_PROJECTS_ID and MC_API_KEY

set -e
BUILD_SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
export CASM_BUILD_DIR=$(dirname $BUILD_SCRIPTS_DIR)
. $CASM_BUILD_DIR/build_scripts/install-functions.sh
. $CASM_BUILD_DIR/build_scripts/build_versions.sh
detect_os

if [[ "$CASM_OS_NAME" == "osx" ]]; then
    brew install bash-completion curl libmagic
fi

check_var "CASM_CONDA_DIR" "Location to install conda and conda environments" "$HOME/.local/conda"
check_var "CASM_VERSION" "CASM version (used for naming conda env)" "0.3"
check_var "CASM_PYTHON_VERSION" "Python version"
check_var "CASM_ENV_NAME" "CASM conda environment name" "casm_${CASM_VERSION}_py${CASM_PYTHON_VERSION}"
check_var "CASM_TEST_PROJECTS_DIR" "Location to download CASM_test_projects" ""

bash $CASM_BUILD_DIR/build_scripts/conda-devel.sh
