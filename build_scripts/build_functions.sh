# These functions require (set in build_all.sh):
#   CASM_GIT_ID_USER=<github user id> (which repo to pull and build)
#   CASM_DOCKER_ID_USER=<docker hub user id> (where to push docker images)
#   CASM_CONDA_ID_USER=<anaconda cloud user id> (where to push conda packages)
#   CASM_GIT_DIR=/path/to/local/CASMcode (local path to CASMcode repo)

# Get variables
. $CASM_GIT_DIR/build_scripts/build_variables.sh

# function to build and upload package with recipe at conda-recipes/$1/$2 $3=version $4=buildnumber
#   skips package if file '$CASM_GIT_DIR/conda-recipes/$1/$2result_py$CASM_PYTHON_VERSION"_"$3"."$4/DONE' exists
build_conda_package () {
  RECIPE_DIR=$CASM_GIT_DIR/conda-recipes/$1/$2
  RESULT_DIR=$RECIPE_DIR/result_py$CASM_PYTHON_VERSION"_"$3"."$4
  if [ ! -f "$RESULT_DIR/DONE" ]; then
    echo
    echo "!! building: "$1" "$3" feature="$2" build="$4" python="$CASM_PYTHON_VERSION" ..."
    echo

    BUILD_FLAGS="--override-channels -c $CASM_CONDA_CHANNEL "
    BUILD_FLAGS+="-c defaults -c conda-forge -c prisms-center "
    BUILD_FLAGS+="--python $CASM_PYTHON_VERSION "

    UPLOAD_FLAGS="--user $CASM_CONDA_ID_USER "
    UPLOAD_FLAGS+="--label $CASM_CONDA_LABEL "

    mkdir -p $RESULT_DIR \
      && conda build $BUILD_FLAGS $RECIPE_DIR > $RESULT_DIR/tmp.out \
      && LOCATION=$(grep 'conda upload' $RESULT_DIR/tmp.out | cut -f3 -d ' ') \
      && if [ -z $LOCATION ]; then \
           LOCATION=$(grep 'Nothing to test for' $RESULT_DIR/tmp.out | cut -f5 -d ' '); \
         fi \
      && cp $LOCATION $RESULT_DIR \
      && anaconda -t $CASM_CONDA_TOKEN_DIR/conda_api_token upload $UPLOAD_FLAGS $LOCATION \
      && echo "true" > $RESULT_DIR/DONE
  else
    echo "$RESULT_DIR/DONE already exists. skipping..."
  fi
}

# build conda-recipes/*/$1
build_conda_all () {
  ## build each package
  build_conda_package "casm-boost" $1 $CASM_BOOST_VERSION $CASM_BOOST_BUILD_NUMBER \
  && build_conda_package "casm-cpp" $1 $CASM_CONDA_VERSION $CASM_BUILD_NUMBER \
  && build_conda_package "casm" $1 $CASM_CONDA_VERSION $CASM_BUILD_NUMBER
}

# build conda-recipes/casm-python/linux using Docker container $DOCKER_ID_USER/casm-build-condagcc:$CASM_CONDA_VERSION
linux_build_casm_python () {
  CASM_GIT_DIR_INSIDE="/home/casmuser/CASMcode"
  CASM_CONDA_TOKEN_DIR_INSIDE="/home/casmuser/tokens/anaconda"
  docker run --rm -it \
    -e CASM_GIT_ID_USER \
    -e CASM_CONDA_ID_USER \
    -e CASM_GIT_DIR=$CASM_GIT_DIR_INSIDE \
    -e CASM_BUILD_NUMBER=$CASM_BUILD_NUMBER \
    -v $CASM_GIT_DIR:$CASM_GIT_DIR_INSIDE \
    -e CASM_CONDA_TOKEN_DIR=$CASM_CONDA_TOKEN_DIR_INSIDE \
    -v $CASM_CONDA_TOKEN_DIR:$CASM_CONDA_TOKEN_DIR_INSIDE \
    $CASM_DOCKER_ID_USER/casm-build-condagcc:$CASM_BRANCH \
    bash -c "cd $CASM_GIT_DIR_INSIDE && . build_scripts/build_functions.sh && build_conda_package \"casm-python\" \"linux\" $CASM_CONDA_VERSION $CASM_BUILD_NUMBER"
}

# build conda-recipes/*/$1 using Docker container $DOCKER_ID_USER/casm-build-$1:$CASM_CONDA_VERSION
linux_build_conda () {
  CASM_GIT_DIR_INSIDE="/home/casmuser/CASMcode"
  CASM_CONDA_TOKEN_DIR_INSIDE="/home/casmuser/tokens/anaconda"
  docker run --rm -it \
    -e CASM_GIT_ID_USER \
    -e CASM_CONDA_ID_USER \
    -e CASM_GIT_DIR=$CASM_GIT_DIR_INSIDE \
    -v $CASM_GIT_DIR:$CASM_GIT_DIR_INSIDE \
    -e CASM_CONDA_TOKEN_DIR=$CASM_CONDA_TOKEN_DIR_INSIDE \
    -v $CASM_CONDA_TOKEN_DIR:$CASM_CONDA_TOKEN_DIR_INSIDE \
    $CASM_DOCKER_ID_USER/casm-build-$1:$CASM_BRANCH \
    bash -c "cd $CASM_GIT_DIR_INSIDE && . build_scripts/build_functions.sh && build_conda_all \"$1\""
}

# build conda-recipes/casm-python/osx
osx_build_casm_python () {
  export NCPUS=$CASM_OSX_NCPUS
  build_conda_package "casm-python" "osx" $CASM_CONDA_VERSION $CASM_BUILD_NUMBER
  export NCPUS=$CASM_DOCKER_NCPUS
}

# build conda-recipes/*/$1
osx_build_conda () {
  export NCPUS=$CASM_OSX_NCPUS
  build_conda_all $1
  export NCPUS=$CASM_DOCKER_NCPUS
}

# build Docker image $DOCKER_ID_USER/$1:$CASM_BRANCH_$2 from docker/$1/$2/Dockerfile,
# use a third argument to pass other args
build_docker () {
  echo "docker build $CASM_GIT_DIR/docker/$1/$2 $3 -t $CASM_DOCKER_ID_USER/$1:$CASM_BRANCH"_"$2"
  echo "\$3: "$3
  docker build $CASM_GIT_DIR/docker/$1/$2 $3 -t $CASM_DOCKER_ID_USER/$1:$CASM_BRANCH"_"$2
}

# push Docker image from $DOCKER_ID_USER/$1:$CASM_BRANCH"_"$2
push_docker () {
  docker push $CASM_DOCKER_ID_USER/$1:$CASM_BRANCH"_"$2
}

# build and push docker image $DOCKER_ID_USER/$1:$CASM_BRANCH"_"$2 from docker/$1/$2/Dockerfile
build_and_push_docker () {
  build_docker "$1" "$2" "$3"
  push_docker "$1" "$2"
}
