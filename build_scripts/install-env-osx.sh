# install conda packages into CASM_ENV_NAME environment
set -e

check_var "CASM_CONDA_DIR" "Conda location"
check_var "CASM_ENV_NAME" "Conda environment name" "casm"
check_var "CASM_DEPENDENCY_VERSIONS" "(Must have sourced build_scripts/build_versions.sh)"

. $CASM_CONDA_DIR/etc/profile.d/conda.sh

if conda list -n $CASM_ENV_NAME > /dev/null 2>&1; then
  echo "conda environment '$CASM_ENV_NAME' already exists"
  echo "To recreate it, first remove it with: "
  echo "  conda remove --name $CASM_ENV_NAME --all -y"
else
  # install conda packages into CASM_ENV_NAME environment
  conda create -c prisms-center -c bpuchala/label/dev -c defaults -c conda-forge -y -n $CASM_ENV_NAME \
    "python =$CASM_PYTHON_VERSION" \
    "casm-boost =$CASM_BOOST_VERSION $CASM_BOOST_XCODE_BUILD_STR" \
    "m4 >=1.4.18" \
    autoconf \
    automake \
    make \
    libtool \
    pkg-config \
    wget \
    bzip2 \
    six
fi
