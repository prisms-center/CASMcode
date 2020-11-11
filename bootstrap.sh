#!/bin/bash

python3 ./make_Makemodule.py

if [ ! $# -eq 0  ]; then
    echo "Applying specified tag to version string: $1"

    version=$(cat casm_version.txt)
    reversion=$1-${version#*-}
    echo $reversion > casm_version.txt
fi

autoreconf -ivf
