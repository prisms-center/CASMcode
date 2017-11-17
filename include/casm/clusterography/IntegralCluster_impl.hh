#ifndef CASM_IntegralCluster_impl
#define CASM_IntegralCluster_impl

#include "casm/clusterography/IntegralCluster.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"

namespace CASM {

  /// \brief Iterate over all sites in an orbit and insert a UnitCellCoord
  ///
  /// \param orbit an Orbit<IntegralCluster>
  /// \param result an OutputIterator for UnitCellCoord
  ///
  /// \result the resulting OutputIterator
  ///
  /// This simply outputs all UnitCellCoord in all equivalent clusters
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename OutputIterator>
  OutputIterator local_orbit_neighborhood(
    const LocalOrbit<IntegralCluster> &orbit,
    OutputIterator result) {

    for(const auto &equiv : orbit) {
      for(const auto &site : equiv) {
        *result++ = site;
      }
    }
    return result;
  }

  /// \brief Iterate over all sites in all orbits and insert a UnitCellCoord
  ///
  /// \param begin,end Range of Orbit<IntegralCluster>
  /// \param result an OutputIterator for UnitCellCoord
  ///
  /// This simply outputs all UnitCellCoord in all equivalent clusters of each orbit
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename ClusterOrbitIterator, typename OutputIterator>
  OutputIterator local_neighborhood(ClusterOrbitIterator begin, ClusterOrbitIterator end, OutputIterator result) {
    // create a neighborhood of all UnitCellCoord that an Orbitree touches
    for(auto it = begin; it != end; ++it) {
      result = local_orbit_neighborhood(*it, result);
    }
    return result;
  }

  /// \brief Iterate over all sites in an orbit and insert a UnitCellCoord
  ///
  /// \param orbit an Orbit<IntegralCluster>
  /// \param result an OutputIterator for UnitCellCoord
  ///
  /// \result the resulting OutputIterator
  ///
  /// This simply outputs all UnitCellCoord for clusters that include the origin
  /// UnitCell, without any standard order. It uses all clusters that touch origin
  /// unitcell, including translationally equivalent clusters.
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename OutputIterator>
  OutputIterator prim_periodic_orbit_neighborhood(
    const PrimPeriodicOrbit<IntegralCluster> &orbit,
    OutputIterator result) {

    for(const auto &equiv : orbit) {

      // UnitCellCoord for all sites in cluster
      std::vector<UnitCellCoord> coord(equiv.begin(), equiv.end());

      // UnitCellCoord for 'flowertree': all clusters that touch origin unitcell
      //  (includes translationally equivalent clusters)
      for(int ns_i = 0; ns_i < coord.size(); ++ns_i) {
        for(int ns_j = 0; ns_j < coord.size(); ++ns_j) {
          *result++ = UnitCellCoord(coord[ns_j].unit(), coord[ns_j].sublat(), coord[ns_j].unitcell() - coord[ns_i].unitcell());
        }
      }
    }
    return result;
  }

  /// \brief Iterate over all sites in all orbits and insert a UnitCellCoord
  ///
  /// \param begin,end Range of Orbit<IntegralCluster>
  /// \param result an OutputIterator for UnitCellCoord
  ///
  /// This simply outputs all UnitCellCoord for clusters that include the origin
  /// UnitCell, without any standard order. It uses all clusters that touch origin
  /// unitcell, including translationally equivalent clusters.
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename ClusterOrbitIterator, typename OutputIterator>
  OutputIterator prim_periodic_neighborhood(ClusterOrbitIterator begin, ClusterOrbitIterator end, OutputIterator result) {
    // create a neighborhood of all UnitCellCoord that an Orbitree touches
    for(auto it = begin; it != end; ++it) {
      result = prim_periodic_orbit_neighborhood(*it, result);
    }
    return result;
  }

}

#endif
