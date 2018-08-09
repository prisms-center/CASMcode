# source this file to detect bash-completion or fail
detect_os

# check standard locations for bash_completion
if [ "$CASM_OS_NAME" == "osx" ]; then
  if [ -f /usr/local/etc/bash_completion ]; then
    DEFAULT_CASM_BASH_COMPLETION=/usr/local/etc/bash_completion
  fi
elif [ "$CASM_OS_NAME" == "linux" ]; then
  if [ -f /usr/share/bash-completion/bash_completion ]; then
    DEFAULT_CASM_BASH_COMPLETION=/usr/share/bash-completion/bash_completion
  fi
fi

check_var "CASM_BASH_COMPLETION" "bash_completion script location" $DEFAULT_CASM_BASH_COMPLETION
