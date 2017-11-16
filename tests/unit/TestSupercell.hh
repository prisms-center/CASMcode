#ifndef UNIT_TestSupercell
#define UNIT_TestSupercell

#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"

using namespace CASM;

namespace test {

  struct TestSupercell {
    TestSupercell(const PrimClex &primclex, const Eigen::Matrix3i &T);
    TestSupercell(const PrimClex &primclex, const Lattice &lat);

    Supercell scel;
    PrimPeriodicSymCompare<IntegralCluster> prim_sym_compare;
    ScelPeriodicSymCompare<IntegralCluster> scel_sym_compare;
    const SymGroup &scel_fg() const;

    mutable SymGroup m_scel_fg;
  };
}

#endif
