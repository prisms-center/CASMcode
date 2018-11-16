# Install latest miniconda, in user space by default
set -e

detect_os
check_var "CASM_CONDA_DIR" "Location to install conda and conda environments" "$HOME/.local/conda"
check_var "CASM_PYTHON_VERSION" "Default Python version" "3.6"
check_var "CASM_CONDA_BUILD_DIR" "Temporary build directory" "/tmp"

# optionally (if CASM_CONDA_DIST is non-zero length)
#   includes conda-build and anaconda-build for building and uploading new conda packages,
#   by default, does not install
check_var "CASM_CONDA_DIST" "If non-zero length string, install conda-build and anaconda-build" ""

check_program curl

use_conda_msg () {
  echo "# To use conda do: "
  echo "#"
  echo "#     $ . $CASM_CONDA_DIR/etc/profile.d/conda.sh"
  echo "#"
  echo "# To make conda always available add: "
  echo "#"
  echo "#     export CONDA_DIR=\${CONDA_DIR:-\"$CASM_CONDA_DIR\"}"
  echo "#     . \$CONDA_DIR/etc/profile.d/conda.sh"
  echo "#"
  echo "# to your ~/.bashrc or ~/.bash_profile file."
  echo "#"
}

# Check OSTYPE to get correct Miniconda version
if [[ "$CASM_OS_NAME" == "osx" ]]; then
  CASM_MINICONDA_NAME="MacOSX-x86_64"
else
  CASM_MINICONDA_NAME="Linux-x86_64"
fi

if conda_exists; then
  if [[ "$CASM_CONDA_DIR" == "$CASM_FOUND_CONDA_DIR" ]]; then
    echo "#"
    echo "# conda is already installed at: $CASM_CONDA_DIR"
    use_conda_msg
  else
    echo "#"
    echo "# conda is already installed at: " $CASM_FOUND_CONDA_DIR
    echo "#"
    echo "# Set CASM_CONDA_DIR=$CASM_FOUND_CONDA_DIR to confirm using this conda installation."
    echo "# Or, remove '. CASM_FOUND_CONDA_DIR/etc/profile.d/conda.sh' from your ~/.bash_profile file."
    echo "#"
    echo "# Installation stopped."
    echo "#"
    exit 1
  fi
elif [ -f $CASM_CONDA_DIR/bin/conda ]; then
  echo "#"
  echo "# conda is already installed at: $CASM_CONDA_DIR"
  use_conda_msg
else
  curl -sSL https://repo.continuum.io/miniconda/Miniconda${CASM_PYTHON_VERSION:0:1}-latest-$CASM_MINICONDA_NAME.sh -o $CASM_CONDA_BUILD_DIR/miniconda.sh
  mkdir -p $CASM_CONDA_DIR
  bash $CASM_CONDA_BUILD_DIR/miniconda.sh -bfp $CASM_CONDA_DIR
  . $CASM_CONDA_DIR/etc/profile.d/conda.sh
  rm -rf $CASM_CONDA_BUILD_DIR/miniconda.sh
  conda install -y "python =$CASM_PYTHON_VERSION"
  if [ -n "$CASM_CONDA_DIST" ]; then
    conda install -y conda-build conda-verify anaconda-client
  fi
  conda update --all --yes
  conda clean --all --yes

  echo "#"
  echo "# conda has been installed at: $CASM_CONDA_DIR"
  use_conda_msg
fi
