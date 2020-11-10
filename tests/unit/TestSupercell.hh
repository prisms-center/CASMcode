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

  /// AxBxC standard FCC cell
  struct TestStandardFCCSupercell : test::TestSupercell {
    TestStandardFCCSupercell(const PrimClex &fcc_primclex, Index A, Index B, Index C) :
      test::TestSupercell(fcc_primclex, lattice(fcc_primclex, A, B, C)) {}

    static Lattice lattice(const PrimClex &primclex, Index A, Index B, Index C) {
      Eigen::Vector3d a, b, c;
      std::tie(a, b, c) = primclex.prim().lattice().vectors();
      Eigen::Vector3d standard_a = c + b - a;
      Eigen::Vector3d standard_b = a - b + c;
      Eigen::Vector3d standard_c = a + b - c;

      return Lattice((1.0 * A) * standard_a, (1.0 * B) * standard_b, (1.0 * C) * standard_c);
    }

  };
}

#endif
