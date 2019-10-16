#ifndef CASM_ClusterDecl
#define CASM_ClusterDecl

#include "casm/symmetry/OrbitDecl.hh"

namespace CASM {
  namespace xtal {
    class UnitCellCoord;
  }

  using namespace xtal;

  template<typename T> struct traits;

  template<typename ClusterType> class ClusterInvariants;
  template<typename _Base> class GenericCluster;
  template<typename _Base> class ElementWiseSymCluster;
  template<typename _Base> class GenericCoordCluster;

  template<typename CoordType> class CoordCluster;


  /// \brief Cluster of UnitCellCoord
  ///
  /// \ingroup IntegralCluster
  ///
  typedef CoordCluster<UnitCellCoord> IntegralCluster;

  typedef AperiodicOrbit<IntegralCluster> AperiodicIntegralClusterOrbit;
  typedef LocalOrbit<IntegralCluster> LocalIntegralClusterOrbit;
  typedef PrimPeriodicOrbit<IntegralCluster> PrimPeriodicIntegralClusterOrbit;
  typedef ScelPeriodicOrbit<IntegralCluster> ScelPeriodicIntegralClusterOrbit;
  typedef WithinScelOrbit<IntegralCluster> WithinScelIntegralClusterOrbit;

}

#endif
