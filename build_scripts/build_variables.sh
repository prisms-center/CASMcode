# variables used in conda-recipes/<package>/<os>/meta.yaml files:

. $CASM_GIT_DIR/build_scripts/build_versions.sh

# ncpus
export CASM_DOCKER_NCPUS="2"
export CASM_OSX_NCPUS="4"
export NCPUS=$CASM_DOCKER_NCPUS

# get development version tag, uses latest tag (i.e. "v0.3.1") as a reference
conda_dev_version () {
  S=$(git describe --abbrev=6 --dirty --always --tags )
  if [ ${S:0:1} == "v" ]; then
    S=${S:1}
  fi
  S=${S/-/}
  S=${S//-/+}
  echo $S
}

export CASM_CONDA_VERSION=$(conda_dev_version) # choose $(conda_dev_version) or "X.Y.Z"
export CASM_CONDA_LABEL="dev"  # choose "dev" or "main"
export CASM_CONDA_CHANNEL="$CASM_CONDA_ID_USER/label/$CASM_CONDA_LABEL"
