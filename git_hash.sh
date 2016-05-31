#!/bin/bash

_git_hash()
{
    githash=$(git rev-parse --short HEAD)
    changes=$(git status --porcelain --untracked-files=no)

    if [ ! -z "$changes" ]; then
        changes=".changes"
    fi

    echo ".${githash}${changes}"
}


#For releases make release true and it will suppress the git hash
release=false

if [[ "$release" == false ]] && hash git 2> /dev/null ; then
    echo "$(_git_hash)"
fi
