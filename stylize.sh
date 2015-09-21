#!/bin/bash

echo "Stylizing declarations..."
for file in $(git diff --name-only --staged | grep ".hh$"); do
    astyle -n --options=astyle_options $file
done

echo "Stylizing definitions..."
for file in $(git diff --name-only --staged | grep ".cc$"); do
    astyle -n --options=astyle_options $file
done

echo "Stylizing actual programs..."
for file in $(git diff --name-only --staged | grep ".cpp$"); do
    astyle -n --options=astyle_options $file
done

touch .is_stylized.txt
