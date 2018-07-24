# activate conda environment
source activate casm

# bash-completion
. /usr/share/bash-completion/bash_completion

# build and test
CXXFLAGS="-O3 -DNDEBUG -Wno-deprecated-register -Wno-ignored-attributes -Wno-deprecated-declarations";

CONFIGFLAGS="--prefix=$CONDA_PREFIX "
CONFIGFLAGS+="--with-zlib=$CONDA_PREFIX "
CONFIGFLAGS+="--with-bash-completion-dir=${BASH_COMPLETION_DIR:-$BASH_COMPLETION_COMPAT_DIR} "
CONFIGFLAGS+="--with-boost-libdir=$CONDA_PREFIX/lib"

source $TRAVIS_BUILD_DIR/build_scripts/travis-script.sh

# check
ldd $TRAVIS_BUILD_DIR/.libs/ccasm
