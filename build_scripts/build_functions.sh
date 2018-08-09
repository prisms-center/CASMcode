# These functions require:
#   CASM_GIT_ID_USER=<github user id> (which repo to pull and build)
#   CASM_CONDA_ID_USER=<anaconda cloud user id> (where to push conda packages)
#   CASM_BUILD_DIR=/path/to/local/CASMcode (local path to CASMcode repo)

# Get variables
. $CASM_BUILD_DIR/build_scripts/build_variables.sh

# function to build and upload package with recipe at conda-recipes/$1/$2 $3=version $4=buildnumber
#   skips package if file '$CASM_BUILD_DIR/conda-recipes/$1/$2result_py$CASM_PYTHON_VERSION"_"$3"."$4/DONE' exists
build_conda_package () {
  RECIPE_DIR=$CASM_BUILD_DIR/conda-recipes/$1/$2
  RESULT_DIR=$RECIPE_DIR/result_py$CASM_PYTHON_VERSION"_"$3"."$4
  unset LOCATION
  if [ ! -f "$RESULT_DIR/DONE" ]; then
    echo
    echo "!! building: "$1" "$3" feature="$2" build="$4" python="$CASM_PYTHON_VERSION" ..."
    echo

    BUILD_FLAGS="--override-channels -c $CASM_CONDA_CHANNEL "
    BUILD_FLAGS+="-c defaults -c conda-forge -c prisms-center "
    BUILD_FLAGS+="--python $CASM_PYTHON_VERSION "

    UPLOAD_FLAGS="--user $CASM_CONDA_ID_USER "
    UPLOAD_FLAGS+="--label $CASM_CONDA_LABEL "

    echo "!! before conda build !!"
    echo "BUILD_FLAGS: $BUILD_FLAGS"
    echo "UPLOAD_FLAGS: $UPLOAD_FLAGS"
    printenv

    echo "!! begin conda build !!"
    mkdir -p $RESULT_DIR \
      && conda build $BUILD_FLAGS $RECIPE_DIR > $RESULT_DIR/tmp.out \
      || do_if_failed "echo \"!! conda build failed !!\"" \
      || do_if_failed "cat $RESULT_DIR/tmp.out" \
      || print_msg_if_failed "!! conda build failed !!" \
      || return 1

    LOCATION=$(grep 'conda upload' $RESULT_DIR/tmp.out | cut -f3 -d ' ') \
      && if [ -z $LOCATION ]; then \
           LOCATION=$(grep 'Nothing to test for' $RESULT_DIR/tmp.out | cut -f5 -d ' '); \
         fi \
      && cp $LOCATION $RESULT_DIR

    conda create -n casm-anaconda-upload "python =$CASM_PYTHON_VERSION" anaconda-client -y \
      || echo "casm-anaconda-upload already exists"
    conda activate casm-anaconda-upload

    if [ -n "$CASM_CONDA_DRY_RUN" ]; then
      echo "(dry run) anaconda -t $CASM_CONDA_TOKEN upload $UPLOAD_FLAGS $LOCATION"
    else
      (anaconda -t $CASM_CONDA_TOKEN upload $UPLOAD_FLAGS $LOCATION && echo "true" > $RESULT_DIR/DONE) \
        || print_msg_if_failed "anaconda upload failed" \
        || return 1
    fi

    conda deactivate

  else
    echo "$RESULT_DIR/DONE already exists. skipping..."
  fi
}
