#!/bin/bash
# Remove all the generated things not removed by 'make clean'

set -v
set -e

# autotools generated things
find . -name .deps -prune -exec rm -rf {} \;
find . -name *.os -exec rm -f {} \;
find . -name *.dirstamp -exec rm -f {} \;
rm -rf autom4te.cache libtool Makefile Makefile.in aclocal.m4 config.h config.h.in config.log config.status configure stamp-h1
rm -f build-aux/casm_version.txt build-aux/compile build-aux/config.guess build-aux/config.sub build-aux/depcomp build-aux/install-sh build-aux/ltmain.sh build-aux/m4/libtool.m4 build-aux/m4/ltoptions.m4 build-aux/m4/ltsugar.m4 build-aux/m4/ltversion.m4 build-aux/m4/lt~obsolete.m4 build-aux/missing build-aux/test-driver

# casm-python generated things
rm -rf ./python/casm/doc/source/python ./python/casm/doc/build
find . -name __pycache__ -prune -exec rm -rf {} \;

# doxygen genrated things
rm -rf doc/html
