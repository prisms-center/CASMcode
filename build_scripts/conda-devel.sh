# Setup a conda development environment for CASM
set -e
BUILD_SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
export CASM_BUILD_DIR=$(dirname $BUILD_SCRIPTS_DIR)
. $CASM_BUILD_DIR/build_scripts/install-functions.sh

detect_os
check_var "CASM_BUILD_DIR" "CASMcode repository location"
check_var "CASM_CONDA_DIR" "Location to install conda and conda environments" "$HOME/.local/conda"
check_var "CASM_ENV_NAME" "Conda environment name" "casm"

. $CASM_BUILD_DIR/build_scripts/build_versions.sh

# install conda
bash $CASM_BUILD_DIR/build_scripts/install-miniconda.sh \
  || { echo "install miniconda failed"; exit 1; }

# install casm development env
bash $CASM_BUILD_DIR/build_scripts/install-env-$CASM_OS_NAME.sh \
  || { echo "create conda environment '$CASM_ENV_NAME' failed"; exit 1; }

. $CASM_CONDA_DIR/etc/profile.d/conda.sh
conda activate $CASM_ENV_NAME

# create $CONDA_PREFIX/.bash-completion.d and $CONDA_PREFIX/.bash-completion
export CASM_BASH_COMPLETION_DIR=$CONDA_PREFIX/.bash_completion.d
export CASM_WRITE_BC_SHORTCUT=$CONDA_PREFIX/.bash_completion
bash $CASM_BUILD_DIR/build_scripts/install-bash-completion.sh \
  || { echo "setup bash-completion failed"; exit 1; }

# download test projects if able
bash $CASM_BUILD_DIR/build_scripts/install-test-projects.sh \
  || { echo "download test projects failed"; exit 1; }

echo ""
echo "The casm conda development environment has been created."
echo "To activate it do:"
echo "  . $CASM_CONDA_DIR/etc/profile.d/conda.sh"
echo "  conda activate $CASM_ENV_NAME"
