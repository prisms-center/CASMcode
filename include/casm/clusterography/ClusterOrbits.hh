#ifndef CASM_ClusterOrbits
#define CASM_ClusterOrbits

#include <vector>
#include <functional>
#include <utility>
#include <iostream>

#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clusterography/ClusterSymCompare.hh"
#include "casm/clusterography/IntegralClusterSymCompareTraits.hh"


namespace CASM {

  namespace xtal {
    class Site;
    class UnitCellCoord;
  }
  template<typename OrbitType> struct OrbitGenerators;
  class Structure;
  class SymGroup;

  /** \defgroup ClusterOrbits

      \brief Functions and classes for IntegralCluster orbits of any type
      \ingroup Clusterography
      \ingroup IntegralCluster

      "ClusterOrbits" (as opposed to IntegralClusterOrbits) has no dependency on casm/clusterography/ClusterSymCompare.
  */

  typedef std::function<bool (xtal::Site)> SiteFilterFunction;
  typedef std::function<bool (IntegralCluster)> ClusterFilterFunction;
  typedef std::function<std::vector<xtal::UnitCellCoord> (Structure const &, SiteFilterFunction)> CandidateSitesFunction;

  struct IntegralClusterOrbitGenerator {

    IntegralClusterOrbitGenerator(
      const IntegralCluster &_prototype,
      bool _include_subclusters = true);

    IntegralCluster prototype;
    bool include_subclusters;
  };

  /* -- OrbitBranchSpecs Declarations ------------------------------------- */

  /// \brief Store data used to generate an orbit branch of IntegralCluster
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename _OrbitType>
  class OrbitBranchSpecs {

  public:

    typedef std::vector<xtal::UnitCellCoord> Container;
    typedef Container::const_iterator const_iterator;
    typedef _OrbitType OrbitType;
    typedef typename OrbitType::Element ClusterType;
    typedef typename OrbitType::SymCompareType SymCompareType;
    typedef Structure PrimType;

    /// \brief Filter returns true for allowed ClusterType, false for unallowed
    typedef std::function<bool (ClusterType)> Filter;

    /// \brief Constructor
    ///
    /// \param _prim A pointer to const _prim is stored in each cluster
    /// \param _begin,_end Iterators over candidate sites that should be considered
    ///                    when determining prototype clusters
    /// \param _generating_grp Group used by the Orbit constructor
    /// \param _filter Function or functor implementing 'bool filter(ClusterType)', which
    ///                returns false for clusters that should not be used to
    ///                construct an Orbit (i.e. pair distance too large)
    /// \param _sym_compare Functor implementing cluster comparison to facilitate
    ///                     comparison with canonicalization, taking into account symmetry
    ///
    template<typename SiteIterator>
    OrbitBranchSpecs(const PrimType &_prim,
                     SiteIterator _begin,
                     SiteIterator _end,
                     const SymGroup &_generating_grp,
                     Filter _filter,
                     const SymCompareType &_sym_compare) :
      m_prim_ptr(&_prim),
      m_candidate_sites(_begin, _end),
      m_generating_grp(&_generating_grp),
      m_filter(_filter),
      m_sym_compare(_sym_compare) {}

    const PrimType &prim() const {
      return *m_prim_ptr;
    }

    std::pair<const_iterator, const_iterator> candidate_sites() const {
      return std::make_pair(m_candidate_sites.cbegin(), m_candidate_sites.cend());
    }

    const SymGroup &generating_group() const {
      return *m_generating_grp;
    }

    Filter filter() const {
      return m_filter;
    }

    const SymCompareType &sym_compare() const {
      return m_sym_compare;
    }


  private:

    const PrimType *m_prim_ptr;

    Container m_candidate_sites;

    const SymGroup *m_generating_grp;

    Filter m_filter;

    SymCompareType m_sym_compare;

  };

  /// \brief Output the neighborhood of UnitCellCoord within max_radius of any sites in unit cell
  template<typename OutputIterator>
  OutputIterator neighborhood(
    Structure const &unit,
    double max_radius,
    SiteFilterFunction site_filter,
    OutputIterator result,
    double xtal_tol);

  /// \brief Output the neighborhood of sites within cutoff_radius of any sites in the phenomenal
  template<typename OutputIterator>
  OutputIterator neighborhood(
    IntegralCluster const &phenomenal,
    double cutoff_radius,
    SiteFilterFunction site_filter,
    OutputIterator result,
    double xtal_tol);

  /// \brief Iterate over all sites in an orbit and insert a UnitCellCoord
  template<typename OrbitType, typename OutputIterator>
  OutputIterator local_orbit_neighborhood(
    const OrbitType &orbit,
    OutputIterator result);

  /// \brief Iterate over all sites in all orbits and insert a UnitCellCoord
  template<typename ClusterOrbitIterator, typename OutputIterator>
  OutputIterator local_neighborhood(ClusterOrbitIterator begin, ClusterOrbitIterator end, OutputIterator result);

  /// \brief Check if periodic images of sites in an orbit are overlapping in supercells defined by the given superlattice transformation matrix
  template<typename OrbitType>
  bool has_local_neighborhood_overlap(
    std::vector<OrbitType> &local_orbits,
    const Eigen::Matrix3i &transf_mat);

  /// \brief Return superlattice transf. matrices for which has_local_neighborhood_overlap is false
  template<typename OrbitType>
  std::vector<Eigen::Matrix3i> viable_supercells(
    std::vector<OrbitType> &local_orbits,
    const std::vector<Eigen::Matrix3i> &transf_mat);


  /* -- Custom clusters --- */

  /// \brief Given a cluster, generate all subcluster generators
  template<typename OrbitType>
  OrbitGenerators<OrbitType> &insert_subcluster_generators(
    typename OrbitType::Element cluster,
    OrbitGenerators<OrbitType> &generators);


  /* -- Orbit access/usage function declarations ------------------------------------- */

  /// \brief Returns the first range containing orbits of the requested orbit branch in the given range of Orbit
  template<typename OrbitIterator>
  std::pair<OrbitIterator, OrbitIterator> orbit_branch(OrbitIterator begin,
                                                       OrbitIterator end,
                                                       Index size);


  /* -- Cluster Orbit generating function declarations ------------------------------------- */

  /// \brief Generate the asymmetric unit, using OrbitBranchSpecs
  template<typename OrbitType, typename OrbitOutputIterator>
  OrbitOutputIterator make_asymmetric_unit(
    const OrbitBranchSpecs<OrbitType> &specs,
    OrbitOutputIterator result,
    std::ostream &status);

  /// \brief Use orbits of size n to generate orbits of size n+1
  template<typename OrbitType, typename OrbitInputIterator, typename OrbitOutputIterator>
  OrbitOutputIterator make_next_orbitbranch(
    OrbitInputIterator begin,
    OrbitInputIterator end,
    const OrbitBranchSpecs<OrbitType> &specs,
    OrbitOutputIterator result,
    std::ostream &status);

  /// \brief Generate orbits of IntegralCluster using OrbitBranchSpecs
  template<typename OrbitBranchSpecsIterator, typename OrbitOutputIterator>
  OrbitOutputIterator make_orbits(
    OrbitBranchSpecsIterator begin,
    OrbitBranchSpecsIterator end,
    const std::vector<IntegralClusterOrbitGenerator> &custom_generators,
    OrbitOutputIterator result,
    std::ostream &status);


  /* -- SymCompareType-Specific IntegralCluster Orbit functions --- */

  /// \brief Generate the asymmetric unit, including all sites
  template<typename OrbitOutputIterator>
  OrbitOutputIterator make_prim_periodic_asymmetric_unit(
    std::shared_ptr<Structure const> prim_ptr,
    SiteFilterFunction const &site_filter,
    double xtal_tol,
    OrbitOutputIterator result,
    std::ostream &status);

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
  /// \param begin,end Range of PrimPeriodicOrbit<IntegralCluster>
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

  /// \brief Iterate over all sites in an orbit and insert a UnitCellCoord
  ///
  /// \param orbit an PrimPeriodicOrbit<IntegralCluster>
  /// \param result an OutputIterator for UnitCellCoord
  ///
  /// \result the resulting OutputIterator
  ///
  /// This simply outputs all UnitCellCoord for clusters that include the origin
  /// UnitCell, without any standard order. It uses all clusters that touch origin
  /// unitcell, including translationally equivalent clusters. Respects translational
  /// properties of local orbits, so can be used when translational type is unknown.
  ///
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename OutputIterator, typename OrbitType>
  OutputIterator flower_neighborhood(
    OrbitType const &orbit,
    OutputIterator result);


  /// \brief Iterate over all sites in all orbits and insert a UnitCellCoord
  ///
  /// \param begin,end Range of PrimPeriodicOrbit<IntegralCluster>
  /// \param result an OutputIterator for UnitCellCoord
  ///
  /// This simply outputs all UnitCellCoord for clusters that include the origin
  /// UnitCell, without any standard order. It uses all clusters that touch origin
  /// unitcell, including translationally equivalent clusters. Respects translational
  /// properties of local orbits, so can be used when translational type is unknown.
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename ClusterOrbitIterator, typename OutputIterator>
  OutputIterator flower_neighborhood(ClusterOrbitIterator begin, ClusterOrbitIterator end, OutputIterator result);

}

#endif
