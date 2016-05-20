#!/bin/bash

VERDIR=./src/casm/version
GITHASH=$(git rev-parse --short HEAD)
CHANGES=$(git status --porcelain --untracked-files=no)

if [ ! -z "$CHANGES" ]; then
    CHANGES=".changes"
fi

sed -e "s|MY_VERSION|$GITHASH$CHANGES|g" < $VERDIR/version_template.cc > $VERDIR/version.cc
