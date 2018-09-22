
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
