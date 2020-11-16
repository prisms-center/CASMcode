## Build and push docker images and conda packages ##

if [[ "$#" -ne 3 ]]; then
    echo "Wrong number of arguments. Expected: "
    echo "  build_conda.sh <github-id> <docker-id> <conda-id>"
    exit 1
fi

set -e

# where to pull CASMcode repo from (https://github.com/$GIT_ID_USER/CASMcode)
export CASM_GIT_ID_USER=$1

# where to push docker images (https://hub.docker.com/u/$DOCKER_ID_USER)
export CASM_DOCKER_ID_USER=$2

# where to push conda packages (https://anaconda.org/$CONDA_ID_USER)
export CASM_CONDA_ID_USER=$3

# the repo top level directory
export CASM_GIT_DIR=$(git rev-parse --show-toplevel)

. $CASM_GIT_DIR/build_scripts/build_functions.sh

if [[ -z "$CASM_CONDA_TOKEN_DIR" ]]; then
    echo "No CASM_CONDA_TOKEN_DIR variable set"
    echo "For passwordless upload, you should set the variable to point"
    echo "at a directory (that is not in the git repository) and do:"
    echo "  anaconda auth --create --name "$CASM_CONDA_ID_USER"_api_token --scopes 'api' > \$CASM_CONDA_TOKEN_DIR/conda_api_token"
    echo "Configure the token for the appropriate channel with"
    echo "  conda config --add channels https://conda.anaconda.org/t/<token>/$CASM_CONDA_ID_USER/label/$CASM_CONDA_LABEL"
    echo "For help, see: https://docs.anaconda.com/anaconda-cloud/user-guide/tasks/work-with-accounts"
    exit 1
fi



### Check and confirm

echo "Will build:"
echo
echo "from: "$CASM_BOOST_URL
echo "  casm-boost "$CASM_BOOST_VERSION" "$CASM_BOOST_BUILD_STR
echo
echo "from: "$CASM_REPO:$CASM_BRANCH
echo "  casm-cpp "$CASM_CONDA_VERSION
echo "  (metapackage) casm "$CASM_CONDA_VERSION
echo
echo "using: "
echo "  c/c++/fortran:"
echo "    "$CASM_XCODE_BUILD_STR" (xcode"$CASM_XCODE_VERSION")"
echo "    "$CASM_CONDAGCC_BUILD_STR" (g[cc|xx|fortran]_linux-64 "$CASM_CONDAGCC_VERSION"*)"
echo
echo "  python:"
echo "    python"$CASM_PYTHON_VERSION"*"
echo
echo "And upload to: "
echo "  anaconda.org/"$CASM_CONDA_ID_USER"/label/"$CASM_CONDA_LABEL
echo
read -p "Continue? [Yy|Nn]: " -n 1 -r
echo    # (optional) move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
  # handle exits from shell or function but don't exit interactive shell
  [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1
fi

### Build base docker images for building linux conda packages
build_docker "casm-build" "condagcc" \
"--build-arg PYTHON_VERSION=$CASM_PYTHON_VERSION "\
"--build-arg CONDAGCC_VERSION=$CASM_CONDAGCC_VERSION"

### Build osx conda packages
if [[ "$OSTYPE" == "darwin"* ]]; then
  osx_build_casm_python
  osx_build_conda "xcode"
fi

### Build linux conda packages
linux_build_conda "condagcc"
linux_build_conda "devtoolset"

### Build casm docker images
build_docker "casm" "condagcc" \
"--build-arg PYTHON_VERSION=$CASM_PYTHON_VERSION "\
"--build-arg CONDAGCC_VERSION=$CASM_CONDAGCC_VERSION"

build_docker "casm" "devtoolset" \
"--build-arg PYTHON_VERSION=$CASM_PYTHON_VERSION "
