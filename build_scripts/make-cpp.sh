# configure && make, with helpful messages

set -e
detect_os

check_var "CASM_BUILD_DIR" "CASMcode repository location"
check_var "CASM_BOOST_PREFIX" "Boost libraries location"
check_var "CASM_BASH_COMPLETION_DIR" "Location to install CASM bash_completion script"
check_var "CASM_CONFIGFLAGS" "configure options"
check_var "CASM_CXXFLAGS" "Compiler flags" "-O3 -Wall -fPIC --std=c++11 -DNDEBUG -Wno-deprecated-register -Wno-ignored-attributes -Wno-deprecated-declarations"
check_var "CASM_CC" "C compiler" ${CC:-"cc"}
check_var "CASM_CXX" "C++ compiler" ${CXX:-"c++"}
check_var "CASM_PYTHON" "Python interpreter" ${PYTHON:-"python"}
check_var "CASM_NCPU" "Compiler -j option" 2
check_var "CASM_MAKE_OPTIONS" "Options to give 'make'" ""

check_program "$CASM_PYTHON"
check_program autoreconf
check_program git
check_program make
check_program "$CASM_CC"
check_program "$CASM_CXX"


### CASM C++ build and test #######################
cd $CASM_BUILD_DIR

# use ccache if available
if which ccache >/dev/null 2>&1; then
  CASM_CC="ccache $CASM_CC"
  CASM_CXX="ccache $CASM_CXX"
fi

echo "version: $(git describe --abbrev=6 --dirty --always --tags)"

echo "git status:"
git status

echo "git diff -–name-only:"
git diff -–name-only

echo "git --no-pager diff:"
git --no-pager diff

echo "git config --global core.autocrlf:"
git config --global core.autocrlf

# C++ and CLI
if ! [ -f ./configure ]; then
  $CASM_PYTHON make_Makemodule.py \
    && ./bootstrap.sh \
    && ./configure CXXFLAGS="${CASM_CXXFLAGS}" CC="$CASM_CC" CXX="$CASM_CXX" PYTHON="$CASM_PYTHON" ${CASM_CONFIGFLAGS} \
    || { echo "configuration failed"; exit 1; }
else
  echo "'configure' already exists; continuing... (delete it to re-configure)"
fi

if grep dirty build-aux/casm_version.txt; then

  echo "version: $(git describe --abbrev=6 --dirty --always --tags)"

  echo "git status:"
  git status

  echo "git diff -–name-only:"
  git diff -–name-only

  echo "git --no-pager diff:"
  git --no-pager diff

  echo "git config --global core.autocrlf:"
  git config --global core.autocrlf

  echo "git is dirty"
  exit 1
else
  echo "git is clean"
fi

echo "make ${CASM_MAKE_OPTIONS:+$CASM_MAKE_OPTIONS} -j $CASM_NCPU"
make ${CASM_MAKE_OPTIONS:+$CASM_MAKE_OPTIONS} -j $CASM_NCPU \
  || { echo "'make ${CASM_MAKE_OPTIONS:+$CASM_MAKE_OPTIONS} -j $CASM_NCPU' failed"; exit 1; }
