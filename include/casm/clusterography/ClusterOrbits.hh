#ifndef CASM_ClusterOrbits
#define CASM_ClusterOrbits

#include "casm/symmetry/Orbit.hh"
#include "casm/clusterography/UnitCellCoordCluster.hh"

namespace CASM {

  /** \defgroup Clusterography

      \brief Functions and classes related to clusters
  */

  /* -- OrbitBranchSpecs Declarations ------------------------------------- */

  /// \brief Store data used to generate an orbit branch of UnitCellCoordCluster
  ///
  ///
  /// \ingroup Clusterography
  ///
  class OrbitBranchSpecs {

  public:

    typedef std::vector<Site> Container;
    typedef Container::const_iterator const_iterator;
    typedef UnitCellCoordCluster::PrimType PrimType;

    /// \brief Filter returns true for allowed UnitCellCoordCluster, false for unallowed
    typedef std::function<bool (UnitCellCoordCluster)> Filter;

    /// \brief Constructor
    ///
    /// \param _prim A pointer to const _prim is stored in each cluster
    /// \param _begin,_end Iterators over candidate sites that should be considered
    ///                    when determining prototype clusters
    /// \param _generating_grp Group used by the Orbit constructor
    /// \param _filter Function or functor implementing 'bool filter(UnitCellCoordCluster)', which
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
                     const SymCompare<ClusterType> &_sym_compare) :
      m_prim_ptr(&_prim),
      m_candidate_sites(_begin, _end),
      m_generating_grp(_generating_grp),
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

    const SymCompare<UnitCellCoordCluster> &sym_compare() const {
      return *m_sym_compare;
    }


  private:

    const PrimType *m_prim_ptr;

    Container m_candidate_sites;

    const SymGroup *m_generating_grp;

    Filter m_filter;

    notstd::cloneable_ptr<SymCompare<UnitCellCoordCluster> > m_sym_compare;

  };


  /* -- Cluster Orbit generating function declarations ------------------------------------- */

  /// \brief Generate the asymmetric unit, including all sites
  template<typename OrbitOutputIterator, typename SymOpIterator>
  OrbitOutputIterator asymmetric_unit(
    const UnitCellCoordCluster::PrimType &prim,
    const SymGroup &generating_grp,
    const SymCompare<UnitCellCoordCluster> &sym_compare,
    OrbitOutputIterator result);

  /// \brief Generate the asymmetric unit, using OrbitBranchSpecs
  template<typename OrbitOutputIterator>
  OrbitOutputIterator asymmetric_unit(const OrbitBrancSpecs &specs, OrbitOutputIterator result);

  /// \brief Use orbits of size n to generate orbits of size n+1
  template<typename ClusterType, typename OrbitInputIterator, typename OrbitOutputIterator>
  OrbitOutputIterator next_orbitbranch(OrbitInputIterator begin,
                                       OrbitInputIterator end,
                                       const OrbitBranchSpecs &specs,
                                       OrbitOutputIterator result,
                                       std::ostream &status);

  /// \brief Generate Orbit<UnitCellCoordCluster>
  template<typename OrbitBranchSpecsIterator, typename OrbitOutputIterator>
  OrbitOutputIterator orbits(OrbitBranchSpecsIterator begin,
                             OrbitBranchSpecsIterator end,
                             OrbitOutputIterator result,
                             std::ostream &status);

  /// \brief Generate Orbit<UnitCellCoordCluster> by specifying max cluster length for each branch
  template<typename OrbitOutputIterator>
  OrbitOutputIterator orbitree(
    const UnitCellCoordCluster::PrimType &prim,
    const SymGroup &generating_grp,
    const std::vector<double> &max_length,
    const std::function<bool (Site)> &site_filter,
    const SymCompare<UnitCellCoordCluster> &sym_compare,
    OrbitOutputIterator result,
    std::stream &status);





  /* -- Orbit access/usage function declarations ------------------------------------- */

  /// \brief Returns the first range containing orbits of the requested orbit branch in the given range of Orbit
  template<typename OrbitIterator>
  std::pair<OrbitIterator, OrbitIterator> orbit_branch(OrbitIterator begin,
                                                       OrbitIterator end,
                                                       unsigned int size);

}

#endif