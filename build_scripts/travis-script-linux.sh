# activate conda environment
source activate casm

# bash-completion setup
export CASM_BASH_COMPLETION_DIR=$CONDA_PREFIX/.bash_completion.d
mkdir -p $CASM_BASH_COMPLETION_DIR
printf "for bcfile in $CONDA_PREFIX/.bash_completion.d/* ; do\n  . \$bcfile\ndone" > $CONDA_PREFIX/.bash_completion
[ -f /usr/share/bash-completion/bash_completion ] && . /usr/share/bash-completion/bash_completion

# build and test variables
CXXFLAGS="-O3 -DNDEBUG -Wno-deprecated-register -Wno-ignored-attributes -Wno-deprecated-declarations";

CONFIGFLAGS="--prefix=$CONDA_PREFIX "
CONFIGFLAGS+="--with-zlib=$CONDA_PREFIX "
CONFIGFLAGS+="--with-boost-libdir=$CONDA_PREFIX/lib "
CONFIGFLAGS+="--with-bash-completion-dir=$CASM_BASH_COMPLETION_DIR "

CASM_NCPU=2

echo "CASM_BASH_COMPLETION_DIR: $CASM_BASH_COMPLETION_DIR"
echo "CXXFLAGS: $CXXFLAGS"
echo "CONFIGFLAGS: $CONFIGFLAGS"

# begin ### build and test #######################
INIT_DIR=$(pwd)
cd $TRAVIS_BUILD_DIR

# C++ and CLI
./bootstrap.sh
./configure CXXFLAGS="${CXXFLAGS}" CC="ccache $CC" CXX="ccache $CXX" ${CONFIGFLAGS}
make -j $CASM_NCPU

make check -j $CASM_NCPU CASM_BOOST_PREFIX="$CONDA_PREFIX"
if [ $TRAVIS_BUILD_DIR/test-suite.log $1 ]; then
  echo "$TRAVIS_BUILD_DIR/test-suite.log: "
  cat $TRAVIS_BUILD_DIR/test-suite.log
else
  echo "does not exist: $TRAVIS_BUILD_DIR/test-suite.log"
fi

# used in casm-python tests
make install

# Python
cd python/casm
pip install .

# skip tests that require a project
unset CASM_TEST_PROJECTS_DIR

# run tests
pytest -r ap -s test_casm

cd $INIT_DIR

# end ### build and test #######################

# source bash-completion scripts
source $CONDA_PREFIX/.bash_completion
