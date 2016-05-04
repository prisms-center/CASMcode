#ifndef CASM_OrbitBranchSpecs
#define CASM_OrbitBranchSpecs

namespace casm {
  
  /** \defgroup Clusterography
      
      \brief Functions and classes related to clusters
  */
  
  /* -- OrbitBranchSpecs Declarations ------------------------------------- */
  
  /// \brief Store data used to generate an orbit branch
  ///
  ///
  /// \ingroup Clusterography
  ///
  template<typename ClusterType>
  class OrbitBranchSpecs {
  
    public:
    
    typedef std::vector<SiteType> Container;
    typedef Container::const_iterator const_iterator; 
    typedef Structure PrimType;
    
    /// \brief Filter returns true for allowed UnitCellCoordCluster, false for unallowed
    typedef std::function<bool (UnitCellCoordCluster)> Filter;
    
    /// \brief Constructor
    ///
    /// \param _prim A pointer to const _prim is stored in each cluster
    /// \param _begin,_end Iterators over candidate sites that should be considered
    ///                    when determining prototype clusters
    /// \param _generating_grp Group passed to the Orbit constructor
    /// \param _filter Function or functor implementing 'bool filter(UnitCellCoordCluster)', which 
    ///                returns false for prototype clusters that should not be used to 
    ///                construct an Orbit (i.e. pair distance too large)
    /// \param _sym_compare Functor implementing cluster comparison to facilitate
    ///                     comparison with canonicalization, taking into account symmetry
    ///
    template<typename SiteIterator>
    OrbitBranchSpecs(const PrimType &_prim,
                     SiteIterator _begin,
                     SiteIterator _end,
                     const Group &_generating_grp,
                     Filter _filter,
                     const SymCompare<ClusterType> &_sym_compare) :
      m_prim_ptr(&_prim),
      m_candidate_sites(_begin, _end),
      m_generating_grp(_generating_grp),
      m_filter(_filter),
      m_sym_compare(_sym_compare) {}
    
    const PrimType& prim() const {
      return *m_prim_ptr;
    }
    
    std::pair<const_iterator, const_iterator> candidate_sites() const {
      return std::make_pair(m_candidate_sites.cbegin(), m_candidate_sites.cend());
    }
    
    const Group &generating_group() const {
      return m_generating_grp;
    }
    
    Filter filter() const {
      return m_filter;
    }
    
    const SymCompare<UnitCellCoordCluster>& sym_compare() const {
      return *m_sym_compare;
    }
    
    
    private:
    
    const PrimType* m_prim_ptr;
    
    Container m_candidate_sites;
    
    Group m_generating_grp;
    
    Filter m_filter;
    
    notstd::cloneable_ptr<SymCompare<UnitCellCoordCluster> > m_sym_compare;
    
  };
  
}

#endif