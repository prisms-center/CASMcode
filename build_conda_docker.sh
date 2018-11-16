# Build and upload conda packages 'casm-cpp' 'casm-python' and 'casm'
#
#   requires variables:
#     CASM_CONDA_TOKEN
#     CASM_CONDA_ID_USER
#     CASM_GIT_ID_USER
#     CASM_BRANCH
#     CASM_DOCKER_CONTAINER
#     CASM_DOCKER_CMD
#     CASM_CONDA_FEATURE (typically xcode/condagcc/condagcc_centos6)
#
#   notable optional env variable:
#     CASM_CONDA_LABEL (default="dev")
#     CASM_BUILD_BOOST (If non-zero length, build casm-boost, otherwise skip)

set -e
export CASM_BUILD_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
. $CASM_BUILD_DIR/build_scripts/install-functions.sh
detect_os

check_var "CASM_DOCKER_CONTAINER" "Docker container used to build the conda packages"
check_var "CASM_CONDA_FEATURE" "Conda feature name (typically xcode/condagcc/condagcc_centos6)"
check_var "CASM_CONDA_LABEL" "Conda channel label (\"dev\" or \"main\")" "dev"
check_var "CASM_BUILD_BOOST" "If non-zero length, build casm-boost, otherwise build casm-python, casm-cpp, and casm" ""
check_var "CASM_GIT_ID_USER" "Pulls CASMcode from this Github user"
check_var "CASM_BRANCH" "Which branch to build"


CASM_BUILD_DIR_INSIDE="/CASMcode"

# ex: CASM_DOCKER_CMD="yum install curl git -y && bash /CASMcode/build_conda.sh"
check_var "CASM_DOCKER_CMD" "Command to run inside the container"

docker run --rm -it \
    -e CASM_GIT_ID_USER \
    -e CASM_CONDA_ID_USER \
    -e CASM_CONDA_TOKEN \
    -e CASM_BRANCH \
    -e CASM_CONDA_FEATURE \
    -e CASM_CONDA_LABEL \
    -e CASM_BUILD_BOOST \
    -v $CASM_BUILD_DIR:$CASM_BUILD_DIR_INSIDE \
    $CASM_DOCKER_CONTAINER /bin/bash -c "$CASM_DOCKER_CMD"
