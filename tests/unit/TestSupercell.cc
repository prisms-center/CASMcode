#include "TestSupercell.hh"
#include "casm/clex/Supercell.hh"

namespace test {

  TestSupercell::TestSupercell(const PrimClex &primclex, const Eigen::Matrix3i &T) :
    scel(&primclex, xtal::make_superlattice(primclex.prim().lattice(), T)),
    prim_sym_compare(primclex.shared_prim(), primclex.crystallography_tol()),
    scel_sym_compare(scel.primclex().shared_prim(), scel.transf_mat(), scel.crystallography_tol()),
    within_scel_sym_compare(scel.primclex().shared_prim(), xtal::make_bring_within_f(scel), scel.crystallography_tol()) {}

  TestSupercell::TestSupercell(const PrimClex &primclex, const Lattice &lat) :
    scel(&primclex, lat),
    prim_sym_compare(primclex.shared_prim(), primclex.crystallography_tol()),
    scel_sym_compare(scel.primclex().shared_prim(), scel.transf_mat(), scel.crystallography_tol()),
    within_scel_sym_compare(scel.primclex().shared_prim(), xtal::make_bring_within_f(scel), scel.crystallography_tol()) {}

  const SymGroup &TestSupercell::scel_fg() const {
    if(!m_scel_fg.size()) {
      m_scel_fg = make_sym_group(scel.sym_info().permute_begin(), scel.sym_info().permute_end(), this->scel.sym_info().supercell_lattice());
    }
    return m_scel_fg;
  }
}

