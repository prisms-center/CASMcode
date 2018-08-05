# default compiler, python, and boost versions

check_var "CASM_GIT_ID_USER" "Pulls CASMcode from this Github user" "prisms-center"

# boost
export CASM_BOOST_VERSION="1.66.0"
export CASM_BOOST_URL="https://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.tar.bz2"
export CASM_BOOST_SHA256="5721818253e6a0989583192f96782c4a98eb6204965316df9f5ad75819225ca9"
export CASM_BOOST_BUILD_NUMBER="0"

# casm
export CASM_BRANCH="0.3.X"
export CASM_REPO="https://github.com/"$CASM_GIT_ID_USER"/CASMcode.git"
export CASM_URL=$CASM_REPO
export CASM_BUILD_NUMBER="0"

# python
export CASM_PYTHON_VERSION=${CASM_PYTHON_VERSION:-"3.6"}

# osx xcode (use system compilers)
export CASM_XCODE_VERSION="9.2"
export CASM_BOOST_XCODE_BUILD_STR="xcode_"$CASM_BOOST_BUILD_NUMBER
export CASM_XCODE_BUILD_STR="xcode_"$CASM_BUILD_NUMBER

# linux condagcc (use conda dist compilers)
export CASM_CONDAGCC_VERSION="7"
export CASM_BOOST_CONDAGCC_BUILD_STR="condagcc_"$CASM_BOOST_BUILD_NUMBER
export CASM_CONDAGCC_BUILD_STR="condagcc_"$CASM_BUILD_NUMBER

# linux devtoolset (use system compilers)
export CASM_DEVTOOLSET_VERSION="7"
export CASM_BOOST_DEVTOOLSET_BUILD_STR="devtoolset_"$CASM_BOOST_BUILD_NUMBER
export CASM_DEVTOOLSET_BUILD_STR="devtoolset_"$CASM_BUILD_NUMBER

# Note that this was sourced
export CASM_DEPENDENCY_VERSIONS="yes"
