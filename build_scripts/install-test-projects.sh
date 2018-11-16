# If CASM_TEST_PROJECTS_DIR set, download CASM test projects data from Materials Commons
#
# If downloading, requires that:
# - secret variable CASM_TEST_PROJECTS_ID is provided
# And either:
# - a $HOME/.materialscommons/config.json file exists
# - or secret variable MC_API_KEY exists
#
# Tests that need test projects should only run if CASM_TEST_PROJECTS_DIR exists
set -e

check_var "CASM_TEST_PROJECTS_DIR" "Location to download CASM_test_projects" ""
check_var "CASM_PYTHON" "Python interpreter" ${PYTHON:-"python"}

echo ""
if [ -n "$CASM_TEST_PROJECTS_DIR" ]; then

  if ! [ -n "$CASM_TEST_PROJECTS_ID" ]; then
    echo "Missing required environment variable CASM_TEST_PROJECTS_ID: Materials Commons CASM_test_projects project id"
    echo "Ask the CASM developers for access at casm-developers@lists.engr.ucsb.edu"
    echo "Or, do 'unset CASM_TEST_PROJECTS_ID' to skip downloading test projects"
    echo ""
    exit 1
  fi

  if ! which mc >/dev/null 2>&1; then
    echo "No 'mc' program detected. Will install materials-commons."
    check_program pip
    pip install materials-commons pathlib2 requests python-magic tabulate sortedcontainers
  fi

  if ! [ -f "$HOME/.materialscommons/config.json" ]; then
    if ! [ -n "$MC_API_KEY" ]; then
      echo "No materialscommons config.json file found."
      echo "Also, missing environment variable MC_API_KEY: Materials Commons API key"
      echo "Ask the CASM developers for access at casm-developers@lists.engr.ucsb.edu"
      echo "Or, do 'unset CASM_TEST_PROJECTS_ID' to skip downloading test projects"
      echo ""
      exit 1
    fi

    $CASM_PYTHON $CASM_BUILD_DIR/build_scripts/write_mc_config.py
  fi

  if [ -d $CASM_TEST_PROJECTS_DIR ]; then
    echo "CASM test projects directory already exists: $CASM_TEST_PROJECTS_DIR"
    echo "To re-download it, first remove it with: "
    echo "  rm -r $CASM_TEST_PROJECTS_DIR"
  else
    echo "Will download CASM test projects"
    (cd $(dirname $CASM_TEST_PROJECTS_DIR) && mc clone $CASM_TEST_PROJECTS_ID) \
      || print_msg_if_failed "Clone test projects failed"

    (cd $CASM_TEST_PROJECTS_DIR && mc down -r 0.3.X) \
      || print_msg_if_failed "Download test project files failed"
  fi
else
  echo "Will not download CASM test projects: CASM_TEST_PROJECTS_DIR not set"
fi
