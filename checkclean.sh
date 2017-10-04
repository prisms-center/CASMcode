#!/bin/bash
rm -f casm_unit_*
find . -name "*.trs" -exec rm {} \;
find . -name 'run_test_*' -not -name '*.in' -exec rm {} \;
rm -rf tests/unit/test_projects

rm -f tests/unit/system/runtime_lib.cc
