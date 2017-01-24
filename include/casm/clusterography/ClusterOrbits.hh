#ifndef CASM_ClusterOrbits
#define CASM_ClusterOrbits

#include <vector>
#include <functional>
#include <utility>
#include <iostream>

#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/clusterography/IntegralCluster.hh"

namespace CASM {

  /*class SymGroup;
  class Index;
  class Structure;
  class Site;*/

  /** \defgroup IntegralCluster

      \brief Functions and classes related to IntegralCluster
      \ingroup Clusterography
      \ingroup CoordCluster
  */

  /* -- OrbitBranchSpecs Declarations ------------------------------------- */

  /// \brief Generate clusters using all Site
  bool all_sites_filter(const Site &site);

  /// \brief Generate clusters using Site with site_occupant.size() > 1
  bool alloy_sites_filter(const Site &site);


  /// \brief Store data used to generate an orbit branch of IntegralCluster
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename _OrbitType>
  class OrbitBranchSpecs {

  public:

    typedef std::vector<UnitCellCoord> Container;
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
  template<typename CoordType, typename OutputIterator>
  OutputIterator neighborhood(
    const Structure &unit,
    double max_radius,
    std::function<bool (CoordType)> site_filter,
    OutputIterator result,
    double xtal_tol);


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

  /// \brief Generate Orbit<IntegralCluster> using OrbitBranchSpecs
  template<typename OrbitBranchSpecsIterator, typename OrbitOutputIterator>
  OrbitOutputIterator make_orbits(
    OrbitBranchSpecsIterator begin,
    OrbitBranchSpecsIterator end,
    const std::vector<IntegralCluster> &custom_generators,
    OrbitOutputIterator result,
    std::ostream &status);


  /* -- Generate prim periodic orbits --------------------------------------- */

  /// \brief Generate the asymmetric unit, including all sites
  template<typename OrbitOutputIterator>
  OrbitOutputIterator make_prim_periodic_asymmetric_unit(
    const Structure &prim,
    const std::function<bool (Site)> &site_filter,
    double xtal_tol,
    OrbitOutputIterator result,
    std::ostream &status);

  /// \brief Generate Orbit<IntegralCluster> by specifying max cluster length for each branch
  template<typename OrbitOutputIterator>
  OrbitOutputIterator make_prim_periodic_orbits(
    const Structure &prim,
    const std::vector<double> &max_length,
    const std::vector<IntegralCluster> &custom_generators,
    const std::function<bool (Site)> &site_filter,
    double xtal_tol,
    OrbitOutputIterator result,
    std::ostream &status);

  /// \brief Generate Orbit<IntegralCluster> from bspecs.json-type JSON input file
  template<typename OrbitOutputIterator>
  OrbitOutputIterator make_prim_periodic_orbits(
    const Structure &prim,
    const jsonParser &bspecs,
    const std::function<bool (Site)> &site_filter,
    double xtal_tol,
    OrbitOutputIterator result,
    std::ostream &status);


  /* -- Orbit access/usage function declarations ------------------------------------- */

  /// \brief Returns the first range containing orbits of the requested orbit branch in the given range of Orbit
  template<typename OrbitIterator>
  std::pair<OrbitIterator, OrbitIterator> orbit_branch(OrbitIterator begin,
                                                       OrbitIterator end,
                                                       Index size);

}

#endif