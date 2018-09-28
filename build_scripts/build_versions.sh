# default compiler, python, and boost versions

set -e

# boost
export CASM_BOOST_VERSION="1.66.0"
export CASM_BOOST_URL="https://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.tar.bz2"
export CASM_BOOST_SHA256="5721818253e6a0989583192f96782c4a98eb6204965316df9f5ad75819225ca9"

# python
export CASM_PYTHON_VERSION=${CASM_PYTHON_VERSION:-"3.6"}

# osx xcode (use system compilers)
export CASM_XCODE_VERSION="9.2"

# linux condagcc (use conda dist compilers)
export CASM_CONDAGCC_VERSION="7"

# linux condagcc (use conda dist compilers)
export CASM_CONDAGCC_CENTOS6_VERSION="7"

# linux devtoolset (use system compilers)
export CASM_DEVTOOLSET_VERSION="7"

# Note that this was sourced
export CASM_DEPENDENCY_VERSIONS="yes"
