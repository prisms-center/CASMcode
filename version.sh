#!/bin/bash

RELEASE="0.0.0"
VERDIR=./src/casm/version
GITHASH=$(git rev-parse --short HEAD)
CHANGES=$(git status --porcelain --untracked-files=no)

if [ ! -z "$CHANGES" ]; then
    CHANGES=".changes"
fi

echo -ne $RELEASE.$GITHASH$CHANGES
#Exchange commented line if release
#echo -ne $RELEASE

#sed -e "s|MY_VERSION|$GITHASH$CHANGES|g" < $VERDIR/version_template.cc > $VERDIR/version.cc
