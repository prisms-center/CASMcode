# variables necessary for travis, linux

check_var "CASM_PREFIX" "Specify the install location"
check_var "CASM_BOOST_PREFIX" "Specify where boost libraries are installed" "$CASM_PREFIX"
check_var "CASM_ZLIB_PREFIX" "Specify where zlib is installed" "$CASM_PREFIX"
check_var "CASM_BASH_COMPLETION" "Specify where bash-completion is installed" "/usr/share/bash-completion/bash_completion"
check_var "CASM_BASH_COMPLETION_DIR" "Specify where to install the casm bash-completion script" "$CASM_PREFIX/.bash_completion.d"

CASM_CONFIGFLAGS="--prefix=$CASM_PREFIX "
CASM_CONFIGFLAGS+="--with-zlib=$CASM_ZLIB_PREFIX "
CASM_CONFIGFLAGS+="--with-boost-libdir=$CASM_BOOST_PREFIX/lib "
CASM_CONFIGFLAGS+="--with-bash-completion-dir=$CASM_BASH_COMPLETION_DIR "
export CASM_CONFIGFLAGS

check_var "CASM_CXXFLAGS" "Compiler flags" "-O3 -Wall -fPIC --std=c++11 -DNDEBUG -Wno-ignored-attributes -Wno-deprecated-declarations -Wno-int-in-bool-context -Wno-sign-compare -Wno-misleading-indentation"
check_var "CASM_CC" "C compiler" ${CC:-"cc"}
check_var "CASM_CXX" "C++ compiler" ${CXX:-"c++"}
check_var "CASM_PYTHON" "Python interpreter" ${PYTHON:-"python"}
check_var "CASM_MAKE_OPTIONS" "Options to give 'make'" ""
