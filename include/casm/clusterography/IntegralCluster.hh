#ifndef CASM_IntegralCluster
#define CASM_IntegralCluster

#include <vector>

#include "casm/clusterography/ClusterDecl.hh"
#include "casm/clusterography/CoordCluster.hh"

namespace CASM {

  /** \defgroup Clusterography

      \brief Functions and classes related to clusters
  */

  /** \defgroup IntegralCluster

      \brief Functions and classes related to IntegralCluster
      \ingroup Clusterography
      \ingroup CoordCluster
  */

  /* -- IntegralCluster Declaration ------------------------------------- */

  /// \brief Print IntegralCluster to stream, using default Printer<IntegralCluster>
  std::ostream &operator<<(std::ostream &sout, const IntegralCluster &clust);

  /// \brief Write IntegralCluster to JSON object
  jsonParser &to_json(const IntegralCluster &clust, jsonParser &json);

  /// \brief Read from JSON
  void from_json(IntegralCluster &clust, const jsonParser &json, double xtal_tol);

  template<>
  struct jsonConstructor<IntegralCluster> {

    /// \brief Construct from JSON
    static IntegralCluster from_json(const jsonParser &json, const Structure &prim, double xtal_tol);
  };



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
    OutputIterator result);

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
  OutputIterator local_neighborhood(ClusterOrbitIterator begin, ClusterOrbitIterator end, OutputIterator result);

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
    OutputIterator result);

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
  OutputIterator prim_periodic_neighborhood(ClusterOrbitIterator begin, ClusterOrbitIterator end, OutputIterator result);

}

#endif
