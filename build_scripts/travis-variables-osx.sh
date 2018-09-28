# variables necessary for travis, osx, building in a conda environment

check_var "CASM_PREFIX" "Specify the install location"
check_var "CASM_BOOST_PREFIX" "Specify where boost libraries are installed" "$CASM_PREFIX"
check_var "CASM_ZLIB_PREFIX" "Specify where zlib is installed" "$CASM_PREFIX"
check_var "CASM_BASH_COMPLETION" "Specify where bash-completion is installed" "/usr/local/etc/bash_completion"
check_var "CASM_BASH_COMPLETION_DIR" "Specify where to install the casm bash-completion script" "$CASM_PREFIX/.bash_completion.d"

CASM_CONFIGFLAGS="--prefix=$CASM_PREFIX "
CASM_CONFIGFLAGS+="--with-zlib=$CASM_ZLIB_PREFIX "
CASM_CONFIGFLAGS+="--with-boost=$CASM_BOOST_PREFIX "
CASM_CONFIGFLAGS+="--with-bash-completion-dir=$CASM_BASH_COMPLETION_DIR "
export CASM_CONFIGFLAGS

check_var "CASM_CXXFLAGS" "Compiler flags" "-O3 -Wall -fPIC --std=c++11 -DNDEBUG -fcolor-diagnostics -Wno-deprecated-register -Wno-ignored-attributes -Wno-deprecated-declarations"
check_var "CASM_CC" "C compiler" ${CC:-"cc"}
check_var "CASM_CXX" "C++ compiler" ${CXX:-"c++"}
check_var "CASM_PYTHON" "Python interpreter" ${PYTHON:-"python"}
check_var "CASM_MAKE_OPTIONS" "Options to give 'make'" ""
