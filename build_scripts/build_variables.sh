# variables used in conda-recipes/<package>/<os>/meta.yaml files:

set -e
. $CASM_BUILD_DIR/build_scripts/build_versions.sh

check_program git

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
export -f conda_dev_version

# choose $(conda_dev_version) or "X.Y.Z"
( cd $CASM_BUILD_DIR; git fetch --unshallow || echo "git fetch --unshallow done"; git tag; )

check_var "CASM_CONDA_VERSION" "Version number for conda package" "$(cd $CASM_BUILD_DIR && conda_dev_version)"

check_var "CASM_CONDA_LABEL" "Conda channel label (\"dev\" or \"main\")" "dev"
check_var "CASM_CONDA_ID_USER" "Where to push conda packages (https://anaconda.org/\$CONDA_ID_USER)"
check_var "CASM_CONDA_CHANNEL" "Conda channel to push package to" "$CASM_CONDA_ID_USER/label/$CASM_CONDA_LABEL"
