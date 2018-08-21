# Linux CASM installation from source, assumes you've installed:
#   g++, boost, git, pip, automake, autoconf, libtool, zlib, bash-completion

if [ "$#" -ne 1 ]; then
    echo "Wrong number of arguments. Expected: "
    echo "  install-casm.sh <prefix>"
    exit 1
fi

set -e  # exit if error

# Variables to set:
CASM_REPO="https://github.com/prisms-center/CASMcode.git"  # where to get CASM
CASM_BRANCH="0.3.X"                                        # which branch/tag to build
CASM_PREFIX=$1                                             # where to install CASM
CASM_BOOST_PREFIX=$CASM_PREFIX  # where to find boost, often same as CASM_PREFIX
NCPUS=4                         # parallel compilation
BUILD_DIR=/tmp                  # where to clone and build

# a local directory to store the casm bash_completion script, use
# your own if you have one
CASM_BASH_COMPLETION_DIR=$CASM_PREFIX/.casm_bash_completion.d

# script to activate user bash_completion scripts,
# after installation, set this environment variable and to activate do:
#   source $CASM_BASH_COMPLETION
CASM_BASH_COMPLETION=$CASM_PREFIX/.casm_bash_completion
mkdir -p $CASM_BASH_COMPLETION_DIR \
  && printf "for bcfile in $CASM_BASH_COMPLETION_DIR/* ; do\n  . \$bcfile\ndone" > $CASM_BASH_COMPLETION

# clone, configure, make, and install CASM
cd $BUILD_DIR \
  && git clone $CASM_REPO \
  && cd CASMcode \
  && git checkout $CASM_BRANCH \
  && pip install --no-cache-dir six \ # necessary for make_Makemodule.py
  && python make_Makemodule.py \ # create Makemodule.am files
  && ./bootstrap.sh \ # create configure script
  && ./configure \
    CXXFLAGS="-O3 -DNDEBUG -Wno-deprecated-register -Wno-ignored-attributes -Wno-deprecated-declarations" \
    --prefix=$CASM_PREFIX \
    --with-boost=$CASM_BOOST_PREFIX \
    #--with-zlib=$CASM_PREFIX \  # necessary when using Conda compiler tools
    #--with-boost-libdir=$CASM_BOOST_PREFIX/lib \  # necessary when using Conda compiler tools
    --with-bash-completion-dir=$CASM_BASH_COMPLETION_DIR \
  && make -j $NCPUS \
  && make install
