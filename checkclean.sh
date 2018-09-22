#!/bin/bash
rm -r ccasm
rm -f casm_unit_*
rm -f .libs/casm_unit_*
find . -name "*.trs" -exec rm {} \;
find . -name 'run_test_*' -not -name '*.in' -exec rm {} \;
find tests -name "*.o" -exec rm {} \;
find tests -name "*.so" -exec rm {} \;
rm -rf tests/unit/test_projects

rm -f tests/unit/system/runtime_lib.cc
