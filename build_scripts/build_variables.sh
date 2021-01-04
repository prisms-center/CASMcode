# variables used in conda-recipes/<package>/<os>/meta.yaml files:

set -e
. $CASM_BUILD_DIR/build_scripts/build_versions.sh

check_program git

# casm
check_var "CASM_GIT_ID_USER" "Pulls CASMcode from this Github user"
check_var "CASM_BRANCH" "Which branch to build"
export CASM_REPO="https://github.com/"$CASM_GIT_ID_USER"/CASMcode.git"
export CASM_URL=$CASM_REPO

# boost
export CASM_BOOST_BUILD_NUMBER="0"

# casm
check_var "CASM_BUILD_NUMBER" "CASM conda build number" "0"

# osx xcode (use system compilers)
export CASM_BOOST_XCODE_BUILD_STR="xcode_"$CASM_BOOST_BUILD_NUMBER
export CASM_XCODE_BUILD_STR="xcode_"$CASM_BUILD_NUMBER

# linux condagcc (use conda dist compilers)
export CASM_BOOST_CONDAGCC_BUILD_STR="condagcc_"$CASM_BOOST_BUILD_NUMBER
export CASM_CONDAGCC_BUILD_STR="condagcc_"$CASM_BUILD_NUMBER

# choose $(conda_version) or "X.Y.Z"
check_var "CASM_CONDA_VERSION" "Version number for conda package"
check_var "CASM_CONDA_LABEL" "Conda channel label (\"dev\" or \"main\")" "dev"
check_var "CASM_CONDA_ID_USER" "Where to push conda packages (https://anaconda.org/\$CONDA_ID_USER)"
check_var "CASM_CONDA_CHANNEL" "Conda channel to push package to" "$CASM_CONDA_ID_USER/label/$CASM_CONDA_LABEL"
