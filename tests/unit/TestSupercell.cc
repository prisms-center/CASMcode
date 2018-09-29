#include "TestSupercell.hh"

namespace test {

  TestSupercell::TestSupercell(const PrimClex &primclex, const Eigen::Matrix3i &T) :
    scel(&primclex, make_supercell(primclex.prim().lattice(), T)),
    prim_sym_compare(primclex),
    scel_sym_compare(scel),
    within_scel_sym_compare(scel) {}

  TestSupercell::TestSupercell(const PrimClex &primclex, const Lattice &lat) :
    scel(&primclex, lat),
    prim_sym_compare(primclex),
    scel_sym_compare(scel),
    within_scel_sym_compare(scel) {}

  const SymGroup &TestSupercell::scel_fg() const {
    if(!m_scel_fg.size()) {
      m_scel_fg = make_sym_group(scel.sym_info().permute_begin(), scel.sym_info().permute_end());
    }
    return m_scel_fg;
  }
}

