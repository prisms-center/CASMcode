#ifndef CASM_ClusterOrbits
#define CASM_ClusterOrbits

namespace casm {
  
  /** \defgroup Clusterography
      
      \brief Functions and classes related to clusters
  */
  
  /* -- Cluster Orbit generating function declarations ------------------------------------- */
  
  /// \brief Return the orbit of empty clusters
  template<ClusterType>
  Orbit<ClusterType> null_orbit(const PrimType &prim);
  
  /// \brief Calculate the asymmetric unit of a PrimMotif
  template<typename OrbitOutputIterator>
  OrbitOutputIterator asymmetric_unit(PrimMotif prim, OrbitOutputIterator result);
  
  /// \brief Calculate the asymmetric unit
  template<typename ClusterType, typename OrbitOutputIterator>
  OrbitOutputIterator asymmetric_unit(const OrbitBranchSpecs<ClusterType> &specs, OrbitOutputIterator result);
  
  /// \brief Use orbits of size n to generate orbits of size n+1
  template<typename ClusterType, typename OrbitInputIterator, typename OrbitOutputIterator>
  OrbitOutputIterator next_orbitbranch(OrbitInputIterator begin,
                                       OrbitInputIterator end,
                                       const OrbitBranchSpecs<ClusterType> &specs,
                                       OrbitOutputIterator result);
  
  /// \brief Get Orbit<UnitCellCoordCluster>
  template<typename OrbitBranchSpecsIterator, typename OrbitOutputIterator>
  OrbitOutputIterator orbits(OrbitBranchSpecsIterator begin,
                             OrbitBranchSpecsIterator end,
                             OrbitOutputIterator result);
  
  
  
  /* -- Orbit access/usage function declarations ------------------------------------- */
  
  /// \brief Returns the first range containing orbits of the requested orbit branch in the given range of Orbit
  template<typename OrbitIterator>
  std::pair<OrbitIterator, OrbitIterator> orbit_branch(OrbitIterator begin, 
                                                       OrbitIterator end, 
                                                       unsigned int size);
  
  
  
  
  
  
  /* -- Cluster Orbit generating function definitions ------------------------------------- */
  
  
  /// \brief Return the orbit of empty clusters
  template<ClusterType>
  Orbit<ClusterType> null_orbit(const PrimType &prim, const SymCompare<ClusterType>& compare) {
    
    ClusterType null(prim);
    
    std::vector<SymOp> ops(1, SymOp());
    SymGroup trivial(ops.cbegin(), ops.cend());
    
    return Orbit<ClusterType>(null, trivial, compare);    
  }
  
  /// \brief Generate the asymmetric unit for a PrimMotif
  ///
  /// \param prim A PrimType
  /// \param result An output iterator for orbits of UnitCellCoordCluster
  ///
  /// \relates UnitCellCoordCluster
  ///
  template<typename OrbitOutputIterator>
  OrbitOutputIterator asymmetric_unit(const PrimType& prim, OrbitOutputIterator result) {}
  
  /// \brief Calculate the asymmetric unit
  ///
  /// \param specs ClusterSpecs
  /// \param result An output iterator for orbits of UnitCellCoordCluster
  ///
  /// \relates UnitCellCoordCluster
  ///
  ///
  template<typename OrbitOutputIterator>
  OrbitOutputIterator asymmetric_unit(const ClusterSpecs &specs, OrbitOutputIterator result) {}
  
  /// \brief Use orbits of size n to generate orbits of size n+1
  ///
  /// \param begin,end A range of input orbits of size n
  /// \param specs ClusterSpecs for orbits of size n+1
  /// \param result An output iterator for orbits of UnitCellCoordCluster
  ///
  /// \relates UnitCellCoordCluster
  ///
  template<typename OrbitInputIterator, typename OrbitOutputIterator>
  OrbitOutputIterator next_orbitbranch(OrbitInputIterator begin,
                                       OrbitInputIterator end,
                                       const ClusterSpecs &specs,
                                       OrbitOutputIterator result) {}
  
  /// \brief Get Orbit<UnitCellCoordCluster>
  ///
  /// \param begin,end Iterators over range of ClusterSpecs, with one for each orbit 
  ///                  branch to be calculated
  /// \param result OutputIterator for resulting Orbit<UnitCellCoordCluster>
  ///
  ///
  /// \relates UnitCellCoordCluster
  ///
  template<typename ClusterSpecsIterator, typename OrbitOutputIterator>
  OrbitOutputIterator orbits(ClusterSpecsIterator begin,
                             ClusterSpecsIterator end,
                             OrbitOutputIterator result) {}
  
  
  /* -- Cluster Orbit access/usage function definitions ------------------------------------- */
  
  
  /// \brief Returns the first range containing orbits of the requested orbit branch in the given range of Orbit
  template<typename OrbitIterator>
  std::pair<OrbitIterator, OrbitIterator> orbit_branch(OrbitIterator begin, 
                                                       OrbitIterator end, 
                                                       unsigned int size) {
    
    auto branch_begin = std::find_if(begin, 
                                     end, 
                                     [&](const typename OrbitIterator::value_type &orbit) {
                                       return orbit.prototype().size() == size;
                                     });
    
    auto branch_end = std::find_if(branch_begin, 
                                   end, 
                                   [&](const typename OrbitIterator::value_type &orbit) {
                                     return orbit.prototype().size() == size+1;
                                   });
    
    return std::pair<OrbitIterator, OrbitIterator>(branch_begin, branch_end);
  }
  
  
  
  
}

#endif