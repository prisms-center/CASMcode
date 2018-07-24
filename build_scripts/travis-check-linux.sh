check_dir () {
  if [ -d $1 ]; then ls $1; else echo "does not exist: $1"; fi
}

check_program () {
  which $1
  if [ $? -eq 0 ]; then
    $1 --version
  else
    echo "does not exist: $1"
  fi
}

check_python () {
  which python
  if [ $? -eq 0 ]; then
    python --version
    python -c "import sys; print(sys.path)"
  else
    echo "does not exist: python"
  fi
}

ldd_check () {
  if [ -f $1 ]; then
    ldd $1
  else
    echo "does not exist: $1"
  if
}

run_checks () {
  printenv
  check_dir $HOME/.ccache
  check_dir $HOME/.local/conda
  check_program conda
  check_program $CXX
  check_program ccache
  check_python
  ldconfig -p
}
