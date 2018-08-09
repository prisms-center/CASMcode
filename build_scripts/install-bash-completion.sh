# bash-completion setup
set -e

check_var "CASM_BASH_COMPLETION_DIR" "Directory where the casm bash-completion script is installed"
check_var "CASM_WRITE_BC_SHORTCUT" "If non-zero length, write bc shortcut script" ""

# Writes a shortcut script at $CASM_WRITE_BC_SHORTCUT so that you can do:
#  `source $CASM_WRITE_BC_SHORTCUT`
# to load all the bash-completion scripts in $CASM_BASH_COMPLETION_DIR
if [ -n "CASM_WRITE_BC_SHORTCUT" ]; then
  printf "for bcfile in $CASM_BASH_COMPLETION_DIR/* ; do\n  . \$bcfile\ndone" > $CASM_WRITE_BC_SHORTCUT
fi

mkdir -p $CASM_BASH_COMPLETION_DIR
