# setup CASM conda build environment
#   if CASM_TEST_PROJECTS_DIR is non-zero length, requires
#     CASM_TEST_PROJECTS_ID and MC_API_KEY

set -e
BUILD_SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
export CASM_BUILD_DIR=$(dirname $BUILD_SCRIPTS_DIR)
export CASM_REPO_SLUG=${CASM_REPO_SLUG:-$TRAVIS_REPO_SLUG}
export CASM_GIT_ID_USER=${CASM_GIT_ID_USER:-${CASM_REPO_SLUG%/*}}
export CASM_BRANCH=${CASM_BRANCH:-$TRAVIS_BRANCH}

### Nothing past here should use travis-ci variables

. $CASM_BUILD_DIR/build_scripts/install-functions.sh
. $CASM_BUILD_DIR/build_scripts/build_versions.sh
detect_os

check_var "CASM_CONDA_DIR" "Location to install conda and conda environments" "$HOME/.local/conda"
check_var "CASM_VERSION" "CASM version (used for naming conda env)" "$CASM_BRANCH"
check_var "CASM_PYTHON_VERSION" "Python version"
check_var "CASM_ENV_NAME" "CASM conda environment name" "casm_${CASM_VERSION}_py${CASM_PYTHON_VERSION}"
check_var "CASM_TEST_PROJECTS_DIR" "Location to download CASM_test_projects" ""

if [[ "$CASM_OS_NAME" == "osx" ]]; then
    brew install bash-completion curl libmagic
fi

bash $CASM_BUILD_DIR/build_scripts/conda-devel.sh
