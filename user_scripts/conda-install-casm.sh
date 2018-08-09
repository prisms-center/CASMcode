### miniconda & CASM auto-install ###

### begin install-functions.sh ###
# Check if variable is of non-zero length, and either use a default value or exit with message
#   $1 variable name as string; $2 description; (optional) $3 default value
check_var () {
  NAME="$1"
  VALUE=${!NAME}
  if ! [ -n "$VALUE" ]; then
    if [[ "$#" -eq 3 ]]; then
      echo "${NAME}=\"$3\" # (default)"
      export "${NAME}"="$3"
    else
      echo "missing required environment variable $1: $2"
      return 1
    fi
  else
    echo "$NAME=\"$VALUE\""
  fi
}
export -f check_var

# Check if variable is of non-zero length, and either use a (non-secret) default value or exit with message
#   $1 variable name as string; $2 description; (optional) $3 default value
check_secret_var () {
  NAME="$1"
  VALUE=${!NAME}
  if ! [ -n "$VALUE" ]; then
    if [[ "$#" -eq 3 ]]; then
      echo "${NAME}=\"$3\" # (default)"
      export "${NAME}"="$3"
    else
      echo "missing required environment variable $1: $2"
      return 1
    fi
  else
    echo "$NAME exists"
  fi
}
export -f check_secret_var

# Require a program be found with "which"
#   $1 program name
check_program () {
  if ! which $1 > /dev/null; then
    echo "missing required program: $1"
    return 1
  fi
}
export -f check_program

# Require a file exists at location
#   $1 location
check_file () {
  if ! [ -f "$1" ]; then
    echo "missing required file: $1"
    return 1
  fi
}
export -f check_file

# Require a directory exists at location
#   $1 location
check_dir () {
  if ! [ -d "$1" ]; then
    echo "missing required directory: $1"
    return 1
  fi
}
export -f check_dir

# Detect OSTYPE and set CASM_OS_NAME to "osx" or "linux"
detect_os () {
  if [[ "$OSTYPE" == "darwin"* ]]; then
    CASM_OS_NAME="osx"
    return 0
  elif [[ "$OSTYPE" == "linux"* ]]; then
    CASM_OS_NAME="linux"
    return 0
  else
    echo "Detected OSTYPE: $OSTYPE"
    echo "CASM is only supported on linux and osx. Exiting..."
    exit 1
  fi
}
export -f detect_os

# Do a command ($1) if previous exit code is not zero, then forward that exit code
#   $1 command
do_if_failed () {
  CODE=$?
  if [ $CODE -ne 0 ]; then
    $1
    return $CODE
  else
    return $CODE
  fi
}
export -f do_if_failed

# echo a string ($1) if previous exit code is not zero, then forward that exit code
#   $1 string to echo
print_msg_if_failed () {
  CODE=$?
  if [ $CODE -ne 0 ]; then
    echo $1
    return $CODE
  else
    return $CODE
  fi
}
export -f print_msg_if_failed

find_conda_dir () {
  unset CASM_FOUND_CONDA_DIR

  if [ -n "$_CONDA_ROOT" ]; then
    CASM_FOUND_CONDA_DIR=$_CONDA_ROOT
  else
    return 1
  fi
}
export -f find_conda_dir

conda_exists () {
  if [[ "$(type -t conda)" == "file" ]]; then
    VERS="$($(which conda) --version 2>&1)"
    LOC="$(dirname $(dirname $(which conda)))"

    echo "# Found: $VERS"
    echo "# Installed here: $LOC"
    echo "#"
    echo "# Your shell has not been properly configured to use 'conda activate'."
    echo "#"
    echo "# If you have not, please update to conda >=4.4."
    echo "# To update conda, do:"
    echo "#"
    echo "#     $ conda update conda"
    echo "#"
    echo "# Previous to conda 4.4, the recommended way to activate conda was to modify PATH in"
    echo "# your ~/.bash_profile file to include '$LOC/bin'. Now to use conda do: "
    echo "#"
    echo "#     $ . $LOC/etc/profile.d/conda.sh"
    echo "#"
    exit 1
  fi
  find_conda_dir
  return $?
}
export -f conda_exists

# typical use: print_source_conda $CASM_CONDA_DIR
print_source_conda () {
  if ! [ -n "$_CONDA_ROOT" ]; then
    echo "#     $ . $1/etc/profile.d/conda.sh"
  else
    exit 1
  fi
}
export -f print_source_conda


### end install-functions.sh ###

### begin install-miniconda.sh ###
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
    conda install -y conda-build anaconda-client
  fi
  conda update --all --yes
  conda clean --all --yes

  echo "#"
  echo "# conda has been installed at: $CASM_CONDA_DIR"
  use_conda_msg
fi
### end install-miniconda.sh ###

set -e

detect_os
check_var "CASM_CONDA_DIR" "Location to install conda and conda environments" "$HOME/.local/conda"
check_var "CASM_PYTHON_VERSION" "Default Python version" "3.6"
check_var "CASM_VERSION" "CASM version to install" "0.3"
check_var "CASM_BUILD_NUMBER" "CASM build number" ""
check_var "CASM_ENV_NAME" "CASM conda environment name" "casm_${CASM_VERSION}_py${CASM_PYTHON_VERSION}"
check_var "CASM_CHANNELS" "Channels for downloading casm & dependencies" "--override-channels -c bpuchala/label/dev -c prisms-center -c defaults -c conda-forge"
check_var "CASM_CONDA_CREATE_FLAGS" "conda create options" "-y"

if ! conda_exists; then
  . $CASM_CONDA_DIR/etc/profile.d/conda.sh
fi

use_casm_msg () {
  echo "# To use CASM do:"
  echo "#"
  echo "#     $ conda activate $CASM_ENV_NAME"
  echo "#"
  echo "# To deactivate the current active environment, do:"
  echo "#"
  echo "#     $ conda deactivate"
  echo "#"
  echo "# To remove CASM, do:"
  echo "#"
  echo "#     $ conda remove --name $CASM_ENV_NAME --all"
  echo "#"
  echo "# To remove CASM and conda, do:"
  echo "#"
  echo "#     $ rm -r $CASM_CONDA_DIR"
  echo "#"
}

if conda list -n $CASM_ENV_NAME > /dev/null 2>&1; then
  echo "#"
  echo "# conda environment '$CASM_ENV_NAME' already exists, it will not be re-installed."
  echo "#"
  echo "# To re-install it, first remove it with: "
  echo "#"
  echo "#     $ conda remove --name $CASM_ENV_NAME --all -y"
  echo "#"
  echo ""
else
  CMD=( conda create -n $CASM_ENV_NAME $CASM_CHANNELS "python =$CASM_PYTHON_VERSION" "casm =$CASM_VERSION $CASM_BUILD_NUMBER" "$CASM_CONDA_CREATE_FLAGS" )
  echo "${CMD[@]}"
  "${CMD[@]}"
fi

echo "#"
echo "# Checking installation..."
conda activate $CASM_ENV_NAME

if ! which casm > /dev/null 2>&1; then
  echo "#"
  echo "# 'casm' not found. Unknown installation failure."
  echo "#"
  exit 1
else
  echo "#"
  echo "# CASM is installed."
  echo "#"
  echo ""
fi

echo "#"
echo "# conda is installed at: $CASM_CONDA_DIR"
use_conda_msg
echo "#"
echo "# CASM is installed in the conda environment: $CASM_ENV_NAME"
use_casm_msg
