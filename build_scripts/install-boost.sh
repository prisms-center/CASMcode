# Build minimal boost libraries necessary for CASM (with default g++ or custom gcc)
# Assumes you've installed:
#   gcc, wget, bzip2

if [ "$#" -ne 1 ]; then
    echo "Wrong number of arguments. Expected: "
    echo "  install-boost.sh <prefix>"
    exit 1
fi

set -e   # exit if error

CASM_BOOST_VERSION=1.66.0       # boost version to install
CASM_BOOST_PREFIX=$1            # where to install boost
NCPUS=4                         # parallel compilation
BUILD_DIR=/tmp                  # where to clone and build
#GXX=/path/to/g++               # if using non-standard gcc set location here

# download, build, and install boost
VERS_NUM=${CASM_BOOST_VERSION//./_} \
  && BOOST_LOC=https://dl.bintray.com/boostorg/release/"$CASM_BOOST_VERSION"/source/boost_"$VERS_NUM".tar.bz2 \
  && cd $BUILD_DIR \
  && wget $BOOST_LOC \
  && tar -xf boost_"$VERS_NUM".tar.bz2 \
  && cd $BUILD_DIR/boost_"$VERS_NUM" \
  && [[ -n "$GXX" ]] && echo "using gcc : : $GXX ;" > tools/build/src/user-config.jam \
  && ./bootstrap.sh \
    --prefix=$PREFIX \
    --with-libraries=system,filesystem,program_options,regex,chrono,timer,test \
  && ./b2 cxxflags="-std=c++11 -O3" -j $NCPUS \
  && ./b2 install -j $NCPUS \
  && rm -rf $BUILD_DIR/boost_"$VERS_NUM".tar.bz2 $BUILD_DIR/boost_"$VERS_NUM"
