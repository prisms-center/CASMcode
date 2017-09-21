#!/bin/bash

stringfile=buid-aux/casm_version.txt
echo "Updating version string in $stringfile"
mkdir -p build-aux
#Be sure to not include a newline!
# git describe --abbrev=6 --dirty --always --tags > build-aux/casm_version.txt
echo "`git symbolic-ref HEAD 2> /dev/null | cut -b 12-`-`git log --pretty=format:\"%h\" -1`" >  build-aux/casm_version.txt
