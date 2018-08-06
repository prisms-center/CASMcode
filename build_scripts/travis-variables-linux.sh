# variables necessary for travis, linux

check_var "CASM_PREFIX" "Specify the install location"

CASM_PREFIX=$CASM_PREFIX
export CASM_BOOST_PREFIX=${CASM_BOOST_PREFIX:-$CASM_PREFIX}

export CASM_BASH_COMPLETION_DIR=${CASM_BASH_COMPLETION_DIR:-$CASM_PREFIX/.bash_completion.d}
export CASM_BASH_COMPLETION=${CASM_BASH_COMPLETION:-"/usr/share/bash-completion/bash_completion"}

CASM_CONFIGFLAGS="--prefix=$CASM_PREFIX "
CASM_CONFIGFLAGS+="--with-zlib=$CASM_PREFIX "
CASM_CONFIGFLAGS+="--with-boost-libdir=$CASM_BOOST_PREFIX/lib "
CASM_CONFIGFLAGS+="--with-bash-completion-dir=$CASM_BASH_COMPLETION_DIR "
export CASM_CONFIGFLAGS

export CASM_CXXFLAGS=${CASM_CXXFLAGS:-"-O3 -Wall -fPIC --std=c++11 -DNDEBUG -Wno-deprecated-register -Wno-ignored-attributes -Wno-deprecated-declarations -Wno-int-in-bool-context -Wno-sign-compare"}
export CASM_CC=${CASM_CC:-"$CC"}
export CASM_CXX=${CASM_CC:-"$CXX"}
