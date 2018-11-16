# Linux install latest miniconda in user space, 
# includes conda-build and anaconda-build for building and uploading new conda packages
# Assumes you've got installed:
#   curl, bzip2

CASM_CONDA_DIR="$HOME/.local/conda"
CASM_PYTHON_VERSION="3"
BUILD_DIR=/tmp

RUN curl -sSL https://repo.continuum.io/miniconda/Miniconda${CASM_PYTHON_VERSION:0:1}-latest-Linux-x86_64.sh -o $BUILD_DIR/miniconda.sh \
  && mkdir -p $CASM_CONDA_DIR \
  && bash $BUILD_DIR/miniconda.sh -bfp $CASM_CONDA_DIR \
  && PATH="$CASM_CONDA_DIR/bin:$PATH" \
  && rm -rf $BUILD_DIR/miniconda.sh \
  && conda install -y \
    "python =$CASM_PYTHON_VERSION" \
    conda-build \
    anaconda-client \
  && conda update --all \
  && conda clean --all --yes
