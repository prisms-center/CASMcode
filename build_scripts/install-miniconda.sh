# Install latest miniconda, in user space by default

# Check OSTYPE to get correct Miniconda version
if [[ "$OSTYPE" == "darwin"* ]]; then
  CASM_MINICONDA_NAME="MacOSX-x86_64"
else
  CASM_MINICONDA_NAME="Linux-x86_64"
fi

# Location to install conda and conda environments, defaults to '$HOME/.local/conda'
CASM_CONDA_DIR="${CASM_CONDA_DIR:-$HOME/.local/conda}"

# Python version to use by default, defaults to '3'
CASM_PYTHON_VERSION="${CASM_PYTHON_VERSION:-3}"

# Temporary build directory, defaults to '/tmp'
CASM_CONDA_BUILD_DIR="${CASM_CONDA_BUILD_DIR:-/tmp}"

# optionally (if CASM_CONDA_DIST is non-zero length)
#   includes conda-build and anaconda-build for building and uploading new conda packages,
#   by default, does not install
CASM_CONDA_DIST="${CASM_CONDA_DIST:-}"

if which conda; then
  echo "conda is already installed and in your path at: " $(which conda)
  exit
elif [ -f $CASM_CONDA_DIR/bin/conda ]; then
  echo "conda is already installed at: $CASM_CONDA_DIR/bin/conda"
  echo "to use conda do:"
  echo "  export PATH=$CASM_CONDA_DIR/bin:\$PATH"
  exit
fi

set -v
curl -sSL https://repo.continuum.io/miniconda/Miniconda${CASM_PYTHON_VERSION:0:1}-latest-$CASM_MINICONDA_NAME.sh -o $CASM_CONDA_BUILD_DIR/miniconda.sh
mkdir -p $CASM_CONDA_DIR
bash $CASM_CONDA_BUILD_DIR/miniconda.sh -bfp $CASM_CONDA_DIR
PATH="$CASM_CONDA_DIR/bin:$PATH"
rm -rf $CASM_CONDA_BUILD_DIR/miniconda.sh
conda install -y "python =$CASM_PYTHON_VERSION"
if [ -n "$CASM_CONDA_DIST" ]; then
  conda install -y conda-build anaconda-client
fi
conda update --all --yes
conda clean --all --yes
set +v

if [ -f $CASM_CONDA_DIR/bin/conda ]; then
  echo "conda has been installed at: $CASM_CONDA_DIR/bin/conda"
  echo "to use conda do:"
  echo "  export PATH=$CASM_CONDA_DIR/bin:\$PATH"
else
  echo "installation failed"
  exit 1
fi
