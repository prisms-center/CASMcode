#ifndef CASM_LocalClusters
#define CASM_LocalClusters

#include <iostream>
#include <iomanip>

#include "casm/clusterography/PeriodicClusters.hh"

namespace casm {

  namespace SiteCluster_impl {
  
    template<typename UnitCellSiteOutputIterator>
    UnitCellSiteOutputIterator _possible_unitcell_sites_inclusive(const PrimMotif &prim,
                                                                  const SiteCluster &phenom,
                                                                  double max_length,
                                                                  UnitCellSiteOutputIterator result);
    
    template<typename UnitCellSiteOutputIterator>
    UnitCellSiteOutputIterator _possible_unitcell_sites_exclusive(const PrimMotif &prim,
                                                                  const SiteCluster &phenom,
                                                                  double max_length,
                                                                  UnitCellSiteOutputIterator result);
    
  }
  
  
  /// \brief Generate ClusterSpecs for local cluster orbits
  ClusterSpecs MaxLengthAperiodicCSPECS_Inclusive(const SiteCluster &phenom,
                                                  const Group &generating_grp,
                                                  double length, 
                                                  double invariants_tol);
  
  /// \brief Generate ClusterSpecs for local cluster orbits
  ClusterSpecs MaxLengthAperiodicCSPECS_Exclusive(const SiteCluster &phenom,
                                                  const Group &generating_grp,
                                                  double length, 
                                                  double invariants_tol);
  
  /// \brief Shortcut for generating local orbits around a phenomenal SiteCluster
  ///
  /// - Generates Orbit of Sitecluster that do include the sites in the phenomenal SiteCluster
  /// - Orbit generation is performed using a Compare function that returns true if SiteCluster
  ///   are exactly identical (not translations).
  ///
  /// Example use:
  /// \code
  /// // Tolerance used for filtering and comparing
  /// double tol = 1e-5;
  ///
  /// // Construct settings for single, pair, and triplet clusters with 
  /// // max cluster site distance 5.0, 5.0 and 4.0, respectively
  /// LocalClusterSpecs_Inclusive lspecs(clust, cluster_group(clust, tol) {5.0, 5.0, 4.0}, tol);
  ///
  /// // Generate local orbits
  /// std::vector<Orbit<SiteCluster> > lorbits;
  /// casm::orbits(lspecs.cbegin(), lspecs.cend(), std::back_inserter(lorbits));
  ///
  /// \endcode
  ///
  /// \ingroup Clusterography
  ///
  class LocalClusterSpecs_Inclusive {
    
    public:
    
    typedef std::vector<ClusterSpecs>::const_iterator const_iterator;
    
    /// \brief Construct using phenomenal SiteCluster and max cluster size for orbit branches 1+ 
    ///
    /// Example use:
    /// \code
    /// std::vector<double> max_length = {5.0, 5.0, 4.0};
    /// LocalClusterSpecs_Inclusive specs(clust, cluster_group(clust, tol), max_length.cbegin(), max_length.cend(), tol);
    /// std::vector<Orbit<SiteCluster> > orbits;
    /// casm::orbits(specs.cbegin(), specs.cend(), std::back_inserter(orbits));
    /// \endcode
    /// 
    template<typename LengthIterator>
    LocalClusterSpecs_Inclusive(const SiteCluster &_cluster,
                                const Group &_generating_grp,
                                LengthIterator _begin,
                                LengthIterator _end,
                                double _invariants_tol) {
      
      // emplace the ClusterSpecs for the n=1+ orbit branches
      for(auto it = _begin; it != _end; ++it) {
        m_specs.push_back(MaxLengthAperiodicCSPECS_Inclusive(_cluster, _generating_grp, *it, _invariants_tol));
      }
    }
    
    /// \brief Construct using phenomenal SiteCluster and max cluster size for orbit branches 1+ 
    ///
    /// Example use:
    /// \code
    /// LocalClusterSpecs_Exclusive specs(clust, cluster_group(clust, tol) {5.0, 5.0, 4.0}, tol);
    /// std::vector<Orbit<SiteCluster> > orbits;
    /// casm::orbits(specs.cbegin(), specs.cend(), std::back_inserter(orbits));
    /// \endcode
    /// 
    LocalClusterSpecs_Inclusive(const SiteCluster &_cluster,
                                const Group &_generating_grp,
                                std::initializer_list<double> _max_length,
                                double _invariants_tol) :
      LocalClusterSpecs_Inclusive(_cluster, _generating_grp, _max_length.begin(), _max_length.end(), _invariants_tol) {}
    
    
    const_iterator cbegin() const {
      return m_specs.cbegin();
    }
    
    const_iterator cend() const {
      return m_specs.cend();
    }
    
    
    private:
    
    std::vector<ClusterSpecs> m_specs;
  
  };
  
  
  /// \brief Shortcut for generating local orbits around a phenomenal SiteCluster
  ///
  /// - Generates Orbit of Sitecluster that do not include the sites in the phenomenal SiteCluster
  /// - Orbit generation is performed using a Compare function that returns true if SiteCluster
  ///   are exactly identical (not translations).
  ///
  /// Example use:
  /// \code
  /// // Tolerance used for filtering and comparing
  /// double tol = 1e-5;
  ///
  /// // Construct settings for single, pair, and triplet clusters with 
  /// // max cluster site distance 5.0, 5.0 and 4.0, respectively
  /// LocalClusterSpecs_Exclusive lspecs(clust, cluster_group(clust, tol) {5.0, 5.0, 4.0}, tol);
  ///
  /// // Generate local orbits
  /// std::vector<Orbit<SiteCluster> > lorbits;
  /// casm::orbits(lspecs.cbegin(), lspecs.cend(), std::back_inserter(lorbits));
  ///
  /// \endcode
  ///
  /// \ingroup Clusterography
  ///
  class LocalClusterSpecs_Exclusive {
    
    public:
    
    typedef std::vector<ClusterSpecs>::const_iterator const_iterator;
    
    /// \brief Construct using phenomenal SiteCluster and max cluster size for orbit branches 1+ 
    ///
    /// Example use:
    /// \code
    /// std::vector<double> max_length = {5.0, 5.0, 4.0};
    /// LocalClusterSpecs_Exclusive specs(clust, cluster_group(clust, tol), max_length.cbegin(), max_length.cend(), tol);
    /// std::vector<Orbit<SiteCluster> > orbits;
    /// casm::orbits(specs.cbegin(), specs.cend(), std::back_inserter(orbits));
    /// \endcode
    /// 
    template<typename LengthIterator>
    LocalClusterSpecs_Exclusive(const SiteCluster &_cluster,
                                const Group &_generating_grp,
                                LengthIterator _begin,
                                LengthIterator _end,
                                double _invariants_tol) {
      
      // emplace the ClusterSpecs for the n=1+ orbit branches
      for(auto it = _begin; it != _end; ++it) {
        m_specs.push_back(MaxLengthAperiodicCSPECS_Exclusive(_cluster, _generating_grp, *it, _invariants_tol));
      }
    }
    
    /// \brief Construct using phenomenal SiteCluster and max cluster size for orbit branches 1+ 
    ///
    /// Example use:
    /// \code
    /// LocalClusterSpecs_Exclusive specs(clust, cluster_group(clust, tol) {5.0, 5.0, 4.0}, tol);
    /// std::vector<Orbit<SiteCluster> > orbits;
    /// casm::orbits(specs.cbegin(), specs.cend(), std::back_inserter(orbits));
    /// \endcode
    /// 
    LocalClusterSpecs_Exclusive(const SiteCluster &_cluster,
                                const Group &_generating_grp,
                                std::initializer_list<double> _max_length,
                                double _invariants_tol) :
      LocalClusterSpecs_Exclusive(_cluster, _generating_grp, _max_length.begin(), _max_length.end(), _invariants_tol) {}
    
    
    const_iterator cbegin() const {
      return m_specs.cbegin();
    }
    
    const_iterator cend() const {
      return m_specs.cend();
    }
    
    
    private:
    
    std::vector<ClusterSpecs> m_specs;
  
  };
  
  
  // ---- Definitions -------------------------------------------- //
  
  namespace SiteCluster_impl {
    
    /// \brief Output UnitCellSites within max_length of all sites in a phenomenal cluster
    ///
    /// Inclusive of sites in the phenomenal cluster
    template<typename UnitCellSiteOutputIterator>
    UnitCellSiteOutputIterator 
    _possible_unitcell_sites_inclusive(const PrimMotif &prim,
                                       const SiteCluster &phenom,
                                       double max_length,
                                       UnitCellSiteOutputIterator result) {
      
      SiteCluster clust = phenom;
      clust.sort();
      
      auto lambda = [&](const UnitCellSite &site) {
        *result++ = site + clust[0].unitcell();
      };
      
      _possible_unitcell_sites(prim, 
                               max_length, 
                               boost::make_function_output_iterator(lambda));
      
      return result;
    }
    
    // \brief Output UnitCellSites within max_length of all sites in a phenomenal cluster
    ///
    /// Exclusive of sites in the phenomenal cluster
    template<typename UnitCellSiteOutputIterator>
    UnitCellSiteOutputIterator
    _possible_unitcell_sites_exclusive(const PrimMotif &prim,
                                       const SiteCluster &phenom,
                                       double max_length,
                                       UnitCellSiteOutputIterator result) {
      
      SiteCluster clust = phenom;
      clust.sort();
      
      // filter out sites in the cluster
      auto lambda = [&](const UnitCellSite &site) {
        if(std::find(clust.cbegin(), clust.cend(), site + clust[0].unitcell()) == clust.cend()) {
          *result++ = site + clust[0].unitcell();
          //std::cout << "-> " << site + clust[0].unitcell() << std::endl;
        }
      };
      
      _possible_unitcell_sites(prim,
                               max_length, 
                               boost::make_function_output_iterator(lambda));
      
      return result;
    }
    
    /// \brief Filter for SiteCluster formed by combining with a phenomenal cluster with max site distance < length
    inline std::function<bool (SiteCluster)> _local_max_length_filter(const SiteCluster &phenom, double length) {
      return [=](const SiteCluster &cluster) {
        
        // get site displacements
        ClusterInvariants invariants(phenom + cluster);
        
        // if test cluster is too big, skip
        if( invariants.displacement().back() > length) {
          return false;
        }
        return true;
      };
    }
    
  }
  
  /// \brief Generate ClusterSpecs for local cluster orbits
  ///
  /// \param phenom The phenomenal cluster around which orbits are generated
  /// \param generating_grp Group used to generate orbits. Should be consistent with phenom.
  /// \param length The maximum length that sites in the local cluster can be from each 
  ///               other or the sites in phenom
  /// \param tol Tolerance used for filtering and comparing
  ///
  /// - Sites in the phenomenal cluster are included in the local orbits
  ///
  /// \relates LocalClusterSpecs_Inclusive
  ///
  inline ClusterSpecs MaxLengthAperiodicCSPECS_Inclusive(const SiteCluster &phenom,
                                                         const Group &generating_grp,
                                                         double length, 
                                                         double invariants_tol) {
    
    // generates a list of candidate UnitCellSites, given the specified max length of clusters
    std::vector<UnitCellSite> unitcell_sites;
    SiteCluster_impl::_possible_unitcell_sites_inclusive(
      phenom.prim_motif(), phenom, length, std::back_inserter(unitcell_sites));
    
    return ClusterSpecs(phenom.prim_motif(),
                        unitcell_sites.cbegin(),
                        unitcell_sites.cend(),
                        generating_grp,
                        SiteCluster_impl::_local_max_length_filter(phenom, length + invariants_tol),
                        ExactCompare<SiteCluster>(),
                        [&](const SiteCluster &cluster) { return ClusterInvariants(phenom + cluster);},
                        invariants_tol);
  }
  
  /// \brief Generate ClusterSpecs for local cluster orbits
  ///
  /// \param phenom The phenomenal cluster around which orbits are generated
  /// \param generating_grp Group used to generate orbits. Should be consistent with phenom.
  /// \param length The maximum length that sites in the local cluster can be from each 
  ///               other or the sites in phenom
  /// \param tol Tolerance used for filtering and comparing
  ///
  /// - Sites in the phenomenal cluster are not included in the local orbits
  ///
  /// \relates LocalClusterSpecs_Exclusive
  ///
  inline ClusterSpecs MaxLengthAperiodicCSPECS_Exclusive(const SiteCluster &phenom,
                                                         const Group &generating_grp,
                                                         double length, 
                                                         double invariants_tol) {
    
    // generates a list of candidate UnitCellSites, given the specified max length of clusters
    std::vector<UnitCellSite> unitcell_sites;
    SiteCluster_impl::_possible_unitcell_sites_exclusive(
      phenom.prim_motif(), phenom, length, std::back_inserter(unitcell_sites));
    
    //for(auto it = unitcell_sites.cbegin(); it != unitcell_sites.cend(); ++it) {
    //  std::cout << "> " << *it << std::endl;
    //}
    
    return ClusterSpecs(phenom.prim_motif(),
                        unitcell_sites.cbegin(),
                        unitcell_sites.cend(),
                        generating_grp,
                        SiteCluster_impl::_local_max_length_filter(phenom, length + invariants_tol),
                        ExactCompare<SiteCluster>(),
                        [&](const SiteCluster &cluster) { return ClusterInvariants(phenom + cluster);},
                        invariants_tol);
  }
  
  
  
  
  
  
  
}

#endif
