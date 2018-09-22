#!/bin/bash

stringfile=build-aux/casm_version.txt
echo "Updating version string in $stringfile"
mkdir -p build-aux
#Be sure to not include a newline!
git describe --abbrev=6 --dirty --always --tags > build-aux/casm_version.txt
