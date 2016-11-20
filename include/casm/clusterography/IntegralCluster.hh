#ifndef CASM_IntegralCluster
#define CASM_IntegralCluster

#include <vector>

#include "casm/symmetry/Orbit.hh"
#include "casm/crystallography/PrimGrid.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/clusterography/ClusterInvariants.hh"
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

  /// \brief Cluster of UnitCellCoord
  ///
  /// \ingroup IntegralCluster
  ///
  typedef CoordCluster<UnitCellCoord> IntegralCluster;

  /// \brief Write IntegralCluster to JSON object
  jsonParser &to_json(const IntegralCluster &clust, jsonParser &json);

  /// \brief Read from JSON
  void from_json(IntegralCluster &clust, const jsonParser &json, double xtal_tol);

  template<>
  struct jsonConstructor<IntegralCluster> {

    /// \brief Construct from JSON
    static IntegralCluster from_json(const jsonParser &json, const Structure &prim, double xtal_tol);
  };



  /* -- IntegralClusterSymCompare Declaration ------------------------------------- */

  template<typename Derived> class IntegralClusterSymCompare;

  namespace CASM_TMP {

    /// \brief Traits class for any IntegralClusterSymCompare derived class
    ///
    /// \ingroup IntegralCluster
    ///
    template<typename Derived>
    struct traits<IntegralClusterSymCompare<Derived> > {
      typedef Derived MostDerived;
      typedef IntegralCluster Element;
      typedef ClusterInvariants<IntegralCluster> InvariantsType;
    };
  }

  /// \brief CRTP Base class for IntegralClusterSymCompare
  ///
  /// Implements:
  /// - 'compare_impl', via comparison of cluster size and lexicographical_compare of UnitCellCoord
  ///
  /// Does not implement:
  /// - 'prepare_impl'
  ///
  /// The ClusterSymCompare hierarchy:
  /// - SymCompare
  ///   - ClusterSymCompare
  ///     - IntegralClusterSymCompare (implements 'compare_impl')
  ///       - LocalSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///       - PrimPeriodicSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///       - ScelPeriodicSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename Derived>
  class IntegralClusterSymCompare : public ClusterSymCompare<IntegralClusterSymCompare<Derived> > {

  protected:

    IntegralClusterSymCompare(double tol) :
      ClusterSymCompare<IntegralClusterSymCompare<Derived> >(tol) {}

    /// \brief Orders 'prepared' elements
    ///
    /// - Returns 'true' to indicate A < B
    /// - Equivalence is indicated by \code !compare(A,B) && !compare(B,A) \endcode
    /// - Assumes elements are 'prepared' before being compared
    /// Implementation:
    /// - compare A.size(), B.size()
    /// - lexicographical_compare of element in A and B
    bool compare_impl(const IntegralCluster &A, const IntegralCluster &B) const {
      if(A.size() != B.size()) {
        return A.size() < B.size();
      }
      return lexicographical_compare(A.begin(), A.end(), B.begin(), B.end());
    }
  };

  /* -- LocalSymCompare<IntegralCluster> Declaration ------------------------------------- */

  /// \brief Comparisons of IntegralCluster with aperiodic symmetry
  ///
  /// The ClusterSymCompare hierarchy:
  /// - SymCompare
  ///   - ClusterSymCompare
  ///     - IntegralClusterSymCompare (implements 'compare_impl')
  ///       - LocalSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///       - PrimPeriodicSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///       - ScelPeriodicSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///
  /// \ingroup IntegralCluster
  ///
  template<>
  class LocalSymCompare<IntegralCluster> : public IntegralClusterSymCompare<LocalSymCompare<IntegralCluster> > {

  public:

    /// \brief Constructor
    ///
    /// \param tol Tolerance for invariants_compare of site-to-site distances
    ///
    LocalSymCompare(double tol):
      IntegralClusterSymCompare<LocalSymCompare<IntegralCluster> >(tol) {}

  protected:

    friend class SymCompare<ClusterSymCompare<IntegralClusterSymCompare<LocalSymCompare<IntegralCluster> > > >;

    /// \brief Prepare an element for comparison
    ///
    /// - Sorts UnitCellCoord
    IntegralCluster prepare_impl(IntegralCluster obj) const {
      std::sort(obj.begin(), obj.end());
      return obj;
    }

  };

  typedef LocalSymCompare<IntegralCluster> LocalIntegralClusterSymCompare;



  /* -- PrimPeriodicSymCompare<IntegralCluster> Declaration ------------------------------------- */


  /// \brief Comparisons of IntegralCluster with periodic symmetry of the primitive lattice
  ///
  /// The ClusterSymCompare hierarchy:
  /// - SymCompare
  ///   - ClusterSymCompare
  ///     - IntegralClusterSymCompare (implements 'compare_impl')
  ///       - LocalSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///       - PrimPeriodicSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///       - ScelPeriodicSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///
  /// \ingroup IntegralCluster
  ///
  template<>
  class PrimPeriodicSymCompare<IntegralCluster> : public IntegralClusterSymCompare<PrimPeriodicSymCompare<IntegralCluster> > {

  public:

    /// \brief Constructor
    ///
    /// \param tol Tolerance for invariants_compare of site-to-site distances
    ///
    PrimPeriodicSymCompare(double tol):
      IntegralClusterSymCompare<PrimPeriodicSymCompare<IntegralCluster> >(tol) {}

  protected:

    friend class SymCompare<ClusterSymCompare<IntegralClusterSymCompare<PrimPeriodicSymCompare<IntegralCluster> > > >;

    /// \brief Prepare an element for comparison
    ///
    /// - Sorts UnitCellCoord and translates so that obj[0] is in the origin unit cell
    IntegralCluster prepare_impl(IntegralCluster obj) const {
      if(!obj.size()) {
        return obj;
      }
      std::sort(obj.begin(), obj.end());
      return obj - obj[0].unitcell();
    }

  };

  typedef PrimPeriodicSymCompare<IntegralCluster> PrimPeriodicIntegralClusterSymCompare;


  /* -- ScelPeriodicSymCompare<IntegralCluster> Declaration ------------------------------------- */


  /// \brief Comparisons of IntegralCluster with periodic symmetry of a supercell lattice
  ///
  /// The ClusterSymCompare hierarchy:
  /// - SymCompare
  ///   - ClusterSymCompare
  ///     - IntegralClusterSymCompare (implements 'compare_impl')
  ///       - LocalSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///       - PrimPeriodicSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///       - ScelPeriodicSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///
  /// \ingroup IntegralCluster
  ///
  template<>
  class ScelPeriodicSymCompare<IntegralCluster> : public IntegralClusterSymCompare<ScelPeriodicSymCompare<IntegralCluster> > {

  public:

    /// \brief Constructor
    ///
    /// \param prim_grid A prim_grid reference
    /// \param tol Tolerance for invariants_compare of site-to-site distances
    ///
    ScelPeriodicSymCompare(const PrimGrid &prim_grid, double tol):
      IntegralClusterSymCompare<ScelPeriodicSymCompare<IntegralCluster> >(tol),
      m_prim_grid(prim_grid) {}

  protected:

    friend class SymCompare<ClusterSymCompare<IntegralClusterSymCompare<ScelPeriodicSymCompare<IntegralCluster> > > >;

    /// \brief Prepare an element for comparison
    ///
    /// - Sorts UnitCellCoord and translates so that obj[0] is within the supercell
    IntegralCluster prepare_impl(IntegralCluster obj) const {
      if(!obj.size()) {
        return obj;
      }
      std::sort(obj.begin(), obj.end());
      auto trans = obj[0].unitcell() - m_prim_grid.within(obj[0]).unitcell();
      return obj - trans;
    }

    const PrimGrid &m_prim_grid;

  };

  typedef ScelPeriodicSymCompare<IntegralCluster> ScelPeriodicIntegralClusterSymCompare;

  typedef Orbit<IntegralCluster, LocalSymCompare<IntegralCluster> > LocalIntegralClusterOrbit;
  typedef Orbit<IntegralCluster, PrimPeriodicSymCompare<IntegralCluster> > PrimPeriodicIntegralClusterOrbit;
  typedef Orbit<IntegralCluster, ScelPeriodicSymCompare<IntegralCluster> > ScelPeriodicIntegralClusterOrbit;


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
    const LocalIntegralClusterOrbit &orbit,
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
    const PrimPeriodicIntegralClusterOrbit &orbit,
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