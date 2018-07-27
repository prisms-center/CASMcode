# activate conda environment
source activate casm

# bash-completion setup
export CASM_BASH_COMPLETION_DIR=$CONDA_PREFIX/.bash_completion.d
mkdir -p $CASM_BASH_COMPLETION_DIR
printf "for bcfile in $CONDA_PREFIX/.bash_completion.d/* ; do\n  . \$bcfile\ndone" > $CONDA_PREFIX/.bash_completion
[ -f /usr/local/etc/bash_completion ] && . /usr/local/etc/bash_completion

# build and test variables
CXXFLAGS="-O3 -DNDEBUG -Wno-deprecated-register -Wno-ignored-attributes -Wno-deprecated-declarations";

CONFIGFLAGS="--prefix=$CONDA_PREFIX "
CONFIGFLAGS+="--with-boost=$CONDA_PREFIX "
CONFIGFLAGS+="--with-bash-completion-dir=$CASM_BASH_COMPLETION_DIR "

CASM_NCPU=2

echo "CASM_BASH_COMPLETION_DIR: $CASM_BASH_COMPLETION_DIR"
echo "CXXFLAGS: $CXXFLAGS"
echo "CONFIGFLAGS: $CONFIGFLAGS"
echo "CC: $CC"
echo "CXX: $CXX"

### CASM C++ build and test #######################

# C++ and CLI
./bootstrap.sh \
  && ./configure CXXFLAGS="${CXXFLAGS}" CC="ccache $CC" CXX="ccache $CXX" ${CONFIGFLAGS} \
  && make -j $CASM_NCPU \
  || { echo "make failure"; exit 1;}

# workaround because Mac SIP does not allow propagating DYLD_LIBRARY_PATH
if ! make check -j $CASM_NCPU CASM_BOOST_PREFIX="$CONDA_PREFIX"; then
  shopt -s extglob
  for X in .libs/@(casm|ccasm)*; do
    echo $X
    install_name_tool -add_rpath "$CONDA_PREFIX/lib" $X
  done
fi

# Run tests and print output
if ! make check -j $CASM_NCPU CASM_BOOST_PREFIX="$CONDA_PREFIX"; then
  cat $TRAVIS_BUILD_DIR/test-suite.log
  echo "make check failure"
  exit 1
fi


### CASM Python install and test #######################

# Installed CASM is used in casm-python tests
make install \
  || { echo "make install failure"; exit 1;}

# Python tests
pip install $TRAVIS_BUILD_DIR/python/casm \
  && (cd $TRAVIS_BUILD_DIR/python/casm && pytest -r ap -s test_casm)
