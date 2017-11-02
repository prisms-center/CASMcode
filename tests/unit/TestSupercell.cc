#include "TestSupercell.hh"

namespace test {

  TestSupercell::TestSupercell(const PrimClex &primclex, const Eigen::Matrix3i &T) :
    scel(&primclex, make_supercell(primclex.prim().lattice(), T)),
    prim_sym_compare(primclex),
    scel_sym_compare(scel),
    scel_fg(make_sym_group(scel.permute_begin(), scel.permute_end())) {}

  TestSupercell::TestSupercell(const PrimClex &primclex, const Lattice &lat) :
    scel(&primclex, lat),
    prim_sym_compare(primclex),
    scel_sym_compare(scel),
    scel_fg(make_sym_group(scel.permute_begin(), scel.permute_end())) {}

}

