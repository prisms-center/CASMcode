printenv

echo "Test" | casm init

casm settings --set-compile-options "${CXX:=g++} -O3 -Wall -fPIC --std=c++11 -I../../include"

casm composition --select 0

casm sym

#casm clusters -e
#casm clusters -c basis_sets/bset.default/custom_clusters.json -s

casm bset -u

casm enum --supercells --max 10
casm enum --configs --max 6

casm select --set-on
casm query -k comp corr -o query_results.txt

casm super --configname SCEL4_2_2_1_1_1_0/0 --transf_mat M2
casm super --configname SCEL4_2_2_1_1_1_0/0 --transf_mat M3
casm super --configname SCEL4_2_2_1_1_1_0/1 --transf_mat M3


