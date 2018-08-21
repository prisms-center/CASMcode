# for running "make" for development in a conda environment
# - the conda environment must be activated

### initialization - shouldn't need to touch
set -e
export CASM_BUILD_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
. $CASM_BUILD_DIR/build_scripts/install-functions.sh
detect_os
check_var "CONDA_PREFIX" "Must have the conda environment activated"
export CASM_PREFIX=$CONDA_PREFIX

### end initialization ###

### variables - Control how CASM is built  ###

check_var "CASM_CXXFLAGS" "Compiler flags" ""
check_var "CASM_NCPU" "Compiler -j option" 2

# set OS-dependent variable defaults
#   only CASM_CONFIGFLAGS can't be overridden from this script
. $CASM_BUILD_DIR/build_scripts/travis-variables-$CASM_OS_NAME.sh

### end variables ###


bash $CASM_BUILD_DIR/build_scripts/make-cpp.sh
