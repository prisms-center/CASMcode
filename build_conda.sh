# Build and upload conda packages 'casm-cpp' and 'casm'
#
#   requires variables:
#     CASM_CONDA_TOKEN
#     CASM_CONDA_ID_USER
#     CASM_GIT_ID_USER
#     CASM_BRANCH
#     CASM_CONDA_VERSION (for pre-releases, use X.Y.0-beta.I -> X.Y.devI)
#     CASM_BOOST_XCODE_BUILD_NUMBER (for osx)
#     CASM_XCODE_BUILD_NUMBER (for osx)
#     CASM_BOOST_CONDAGCC_BUILD_NUMBER (for linux)
#     CASM_CONDAGCC_BUILD_NUMBER (for linux)
#
#   notable optional env variable:
#     CASM_CONDA_FEATURE (typically xcode/condagcc/condagcc_centos6)
#     CASM_CONDA_LABEL (default="dev")
#     CASM_BUILD_BOOST (If non-zero length, build casm-boost, otherwise skip)

set -e
export CASM_BUILD_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
. $CASM_BUILD_DIR/build_scripts/install-functions.sh
detect_os

# where to pull CASMcode repo from (https://github.com/$GIT_ID_USER/CASMcode)
check_var "CASM_BRANCH" "Branch to build"
check_var "CASM_CONDA_VERSION" "Version number for conda package"
check_var "CASM_GIT_ID_USER" "Where to pull CASMcode from (https://github.com/\$GIT_ID_USER/CASMcode)"
check_var "CASM_CONDA_ID_USER" "Where to push conda packages (https://anaconda.org/\$CONDA_ID_USER)"
check_secret_var "CASM_CONDA_TOKEN" "Token required to upload conda packages"
check_var "CASM_CONDA_LABEL" "Conda channel label (\"dev\" or \"main\")" "dev"
check_var "CASM_CONDA_DIR" "Location where conda is installed" "$HOME/.local/conda"
check_var "CASM_BUILD_BOOST" "If non-zero length, build casm-boost, otherwise build casm-cpp and casm" ""
check_var "CASM_NCPU" "Compiler -j option" 2

export CASM_CONDA_DIST="yes"
bash $CASM_BUILD_DIR/build_scripts/install-miniconda.sh

. $CASM_CONDA_DIR/etc/profile.d/conda.sh
. $CASM_BUILD_DIR/build_scripts/build_functions.sh

if [ "$CASM_OS_NAME" == "osx" ]; then
  CASM_BOOST_BUILD_STR=$CASM_BOOST_XCODE_BUILD_STR
  CASM_DEFAULT_CONDA_FEATURE="xcode"
elif [ "$CASM_OS_NAME" == "linux" ]; then
  CASM_BOOST_BUILD_STR=$CASM_BOOST_CONDAGCC_BUILD_STR
  CASM_DEFAULT_CONDA_FEATURE="condagcc"
fi
check_var "CASM_CONDA_FEATURE" "Conda feature name (typically xcode/condagcc/condagcc_centos6)" "$CASM_DEFAULT_CONDA_FEATURE"

# only build casm-boost occasionally
if [ -n "$CASM_BUILD_BOOST" ]; then
  build_conda_package "casm-boost-cpp17" $CASM_CONDA_FEATURE $CASM_BOOST_VERSION $CASM_BOOST_BUILD_STR
else
  build_conda_package "casm-cpp" $CASM_CONDA_FEATURE $CASM_CONDA_VERSION $CASM_BUILD_STR
  # build_conda_package "casm" $CASM_CONDA_FEATURE $CASM_CONDA_VERSION $CASM_BUILD_STR
fi
