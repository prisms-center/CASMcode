#ifndef CASM_IntegralClusterSymCompareTraits
#define CASM_IntegralClusterSymCompareTraits

#include "casm/clusterography/ClusterSymCompareDecl.hh"
#include "casm/symmetry/OrbitDecl.hh"

namespace CASM {

namespace xtal {
class UnitCellCoord;
}

class ClusterInvariants;
class IntegralCluster;
class SymOp;
// class WithinScelClusterInvariants;

/** \defgroup Clusterography

    \brief Functions and classes related to clusters
*/

/** \defgroup IntegralClusterOrbits

    \brief Traits and functions for specific types of IntegralCluster orbits
    \ingroup Clusterography
    \ingroup IntegralCluster

*/

/// Traits used for:
/// - PrimPeriodicSymCompare<IntegralCluster>
/// - AperiodicSymCompare<IntegralCluster>
/// - ScelPeriodicSymCompare<IntegralCluster>
template <typename SymCompareType>
struct IntegralClusterSymCompareTraits {
  using Element = IntegralCluster;
  using InvariantsType = ClusterInvariants;

  /// Returns clust[0]
  static xtal::UnitCellCoord position(IntegralCluster const &clust,
                                      SymCompareType const &sym_compare);

  static Element copy_apply(SymOp const &op, IntegralCluster const &clust,
                            SymCompareType const &sym_compare);

  static ClusterInvariants make_invariants(IntegralCluster const &clust,
                                           SymCompareType const &sym_compare);
};

template <>
struct traits<AperiodicSymCompare<IntegralCluster>>
    : public IntegralClusterSymCompareTraits<
          AperiodicSymCompare<IntegralCluster>> {};

template <>
struct traits<PrimPeriodicSymCompare<IntegralCluster>>
    : public IntegralClusterSymCompareTraits<
          PrimPeriodicSymCompare<IntegralCluster>> {};

template <>
struct traits<ScelPeriodicSymCompare<IntegralCluster>>
    : public IntegralClusterSymCompareTraits<
          ScelPeriodicSymCompare<IntegralCluster>> {};

typedef AperiodicOrbit<IntegralCluster> AperiodicIntegralClusterOrbit;
typedef LocalOrbit<IntegralCluster> LocalIntegralClusterOrbit;
typedef PrimPeriodicOrbit<IntegralCluster> PrimPeriodicIntegralClusterOrbit;
typedef ScelPeriodicOrbit<IntegralCluster> ScelPeriodicIntegralClusterOrbit;

// /// Traits used for WithinScelSymCompare<IntegralCluster>
// template <>
// struct traits<WithinScelSymCompare<IntegralCluster>> {
//
//   using Element = IntegralCluster;
//   using InvariantsType = WithinScelClusterInvariants;
//
//   static IntegralCluster bring_within(
//     IntegralCluster clust,
//     WithinScelSymCompare<IntegralCluster> const &sym_compare);
//
//   static Element copy_apply(
//     SymOp const &op,
//     IntegralCluster const &clust,
//     WithinScelSymCompare<IntegralCluster> const &sym_compare);
//
//   static WithinScelClusterInvariants make_invariants(
//     IntegralCluster clust,
//     WithinScelSymCompare<IntegralCluster> const &sym_compare);
// };
//
// typedef WithinScelOrbit<IntegralCluster> WithinScelIntegralClusterOrbit;

}  // namespace CASM

#endif
