# check casm-python, using casm built in $CASM_BUILD_DIR/.libs
set -e
detect_os

check_var "CASM_BUILD_DIR" "CASMcode repository location"

### CASM Python install and test #######################
cd $CASM_BUILD_DIR/python/casm

# make sure we can find the made, but not installed, casm
if [[ "$CASM_OS_NAME" == "osx" ]]; then
  export LIBCASM=$CASM_BUILD_DIR/.libs/libcasm.dylib
elif [[ "$CASM_OS_NAME" == "linux" ]]; then
  export LIBCASM=$CASM_BUILD_DIR/.libs/libcasm.so
else
  exit 1
fi
PATH=$CASM_BUILD_DIR/.libs:$PATH

check_program ccasm
pip install -e .
pip install -r test_requirements.txt
pytest -r ap -s test_casm
