#!/bin/bash

for file in `find src | grep "\.cc" | grep -v "external"`
do
    echo $file
    astyle -n --options=astyle_options $file
done

for file in `find tests | grep "\.cc"`
do
    echo $file
    astyle -n --options=astyle_options $file
done

for file in `find tests | grep "\.cpp"`
do
    echo $file
    astyle -n --options=astyle_options $file
done

for file in `find tests | grep "\.hh"`
do
    echo $file
    astyle -n --options=astyle_options $file
done

for file in `find include | grep "\.hh" | grep -v "external"`
do
    echo $file
    astyle -n --options=astyle_options $file
done

touch .is_stylized.txt
