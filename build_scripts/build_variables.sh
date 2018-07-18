# variables used in conda-recipes/<package>/<os>/meta.yaml files:

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
export CASM_PYTHON_VERSION="3.6"

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

# ncups
export CASM_DOCKER_NCPUS="2"
export CASM_OSX_NCPUS="4"
export NCPUS=$CASM_DOCKER_NCPUS

# get development version tag, uses latest tag (i.e. "v0.3.1") as a reference
conda_dev_version () {
  S=$(git describe --abbrev=6 --dirty --always --tags )
  if [ ${S:0:1} == "v" ]; then
    S=${S:1}
  fi
  S=${S/-/}
  S=${S//-/+}
  echo $S
}

export CASM_CONDA_VERSION="0.3.dev132+g5e2a68" #$(conda_dev_version) # choose $(conda_dev_version) or "X.Y.Z"
export CASM_CONDA_LABEL="dev"  # choose "dev" or "main"
export CASM_CONDA_CHANNEL="$CASM_CONDA_ID_USER/label/$CASM_CONDA_LABEL"
    
