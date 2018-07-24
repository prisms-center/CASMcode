set -x

# common script to be called from travis-script-$TRAVIS_OS_NAME.sh
INIT_DIR=$(pwd)
cd $TRAVIS_BUILD_DIR

CASM_NCPU=2
echo "CASM_NCPU: $CASM_NCPU"

# C++ and CLI
echo "CXXFLAGS: "$CXXFLAGS
echo "CONFIGFLAGS: "$CONFIGFLAGS
./bootstrap.sh
./configure CXXFLAGS="${CXXFLAGS}" ${CONFIGFLAGS}
make -j $CASM_NCPU
make check -j $CASM_NCPU
make install

# Python
cd python/casm
pip install .
pytest -r ap -s test_casm

cd $INIT_DIR
