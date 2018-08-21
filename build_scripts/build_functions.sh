# These functions require:
#   CASM_GIT_ID_USER=<github user id> (which repo to pull and build)
#   CASM_CONDA_ID_USER=<anaconda cloud user id> (where to push conda packages)
#   CASM_BUILD_DIR=/path/to/local/CASMcode (local path to CASMcode repo)

# Get variables
. $CASM_BUILD_DIR/build_scripts/build_variables.sh

# function to build and upload package with recipe at conda-recipes/$1/$2
# $1=package name, $2=feature name, $3=version, $4=buildnumber
build_conda_package () {
  RECIPE_DIR=$CASM_BUILD_DIR/conda-recipes/$1/$2

  conda create -n casm-anaconda-upload "python =$CASM_PYTHON_VERSION" anaconda-client conda-verify -y \
    || echo "casm-anaconda-upload already exists"
  conda activate casm-anaconda-upload

  echo
  echo "!! building: "$1" "$3" feature="$2" build="$4" python="$CASM_PYTHON_VERSION" ..."
  echo

  BUILD_FLAGS="--override-channels -c $CASM_CONDA_CHANNEL "
  BUILD_FLAGS+="-c defaults -c conda-forge -c prisms-center "
  BUILD_FLAGS+="--python $CASM_PYTHON_VERSION "

  UPLOAD_FLAGS="--user $CASM_CONDA_ID_USER --label $CASM_CONDA_LABEL "

  LOCATION=$(conda build $BUILD_FLAGS $RECIPE_DIR --output)
  conda build $BUILD_FLAGS $RECIPE_DIR
  anaconda -t $CASM_CONDA_TOKEN upload $UPLOAD_FLAGS $LOCATION

  conda deactivate
}
