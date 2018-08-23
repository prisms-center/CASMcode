# workaround because Mac SIP does not allow propagating DYLD_LIBRARY_PATH
set -e
detect_os

check_var "CASM_OS_NAME" "Determine with detect_os"

if [[ "$CASM_OS_NAME" == "linux" ]]; then
  exit 1
fi

check_var "CASM_BOOST_PREFIX" "Boost libraries location"
check_var "CASM_NCPU" "Compiler -j option" 2


# workaround because Mac SIP does not allow propagating DYLD_LIBRARY_PATH
shopt -s extglob
shopt -s nullglob
files=(.libs/*(lib)@(casm|ccasm)!(*.la|*.lai) )

unset SET_RPATH
for X in "${files[@]}"; do
  echo "check $X"
  if ! otool -l "$X" | grep LC_RPATH >/dev/null 2>&1; then
    echo "did not find LC_RPATH, setting SET_RPATH=\"true\""
    SET_RPATH="true"
  fi
done

if [ -n "$SET_RPATH" ]; then
  echo ""
  echo "^ That probably failed because osx system integrity protection doesn't allow exporting DYLD_LIBRARY_PATH"
  echo "We'll set rpath manually and try again!"
  echo ""

  # set rpath manually to find boost libraries
  for X in "${files[@]}"; do
    echo "install_name_tool -add_rpath $CASM_BOOST_PREFIX/lib $X"
    install_name_tool -add_rpath "$CASM_BOOST_PREFIX/lib" "$X" \
      || { echo "  already set"; }
  done

  # now re-run checks
  echo "CASM_TESTS: '$CASM_TESTS'"
  echo "CASM_TEST_FLAGS: '$CASM_TEST_FLAGS'"
  echo "CASM_BOOST_PREFIX: '$CASM_BOOST_PREFIX'"

  make check -j $CASM_NCPU CASM_BOOST_PREFIX="$CASM_BOOST_PREFIX" ${CASM_TESTS:+TESTS="$CASM_TESTS"} TEST_FLAGS="$CASM_TEST_FLAGS"


else
  echo "^ That probably failed because a test really did fail"
  exit 1
fi
