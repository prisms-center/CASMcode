#include "casm/clusterography/ClusterOrbits_impl.hh"

namespace CASM {

#define  CLUSTER_ORBITS_INST(ITERATOR,INSERTER,ORBIT,SPECSITERTOR) \
  template INSERTER make_asymmetric_unit<INSERTER,typename ORBIT::SymCompareType>( \
    const IntegralCluster::PrimType &prim, \
    const SymGroup &generating_grp, \
    const typename ORBIT::SymCompareType &sym_compare, \
    INSERTER result); \
  \
  template INSERTER make_asymmetric_unit<ORBIT,INSERTER>( \
    const OrbitBranchSpecs<ORBIT> &specs, \
    INSERTER result, \
    std::ostream &status); \
  \
  template INSERTER make_next_orbitbranch<ORBIT,ITERATOR,INSERTER>( \
    ITERATOR begin, \
    ITERATOR end, \
    const OrbitBranchSpecs<ORBIT> &specs, \
    INSERTER result, \
    std::ostream &status); \
  \
  template INSERTER make_orbits<SPECSITERTOR,INSERTER>( \
    SPECSITERTOR begin, \
    SPECSITERTOR end, \
    INSERTER result, \
    std::ostream &status); \
  \
  template INSERTER make_orbits<INSERTER,typename ORBIT::SymCompareType>( \
    const IntegralCluster::PrimType &prim, \
    const SymGroup &generating_grp, \
    const std::vector<double> &max_length, \
    double crystallography_tol, \
    const std::function<bool (Site)> &site_filter, \
    const typename ORBIT::SymCompareType &sym_compare, \
    INSERTER result, \
    std::ostream &status); \
  \
  template INSERTER make_orbits<INSERTER,typename ORBIT::SymCompareType>( \
    const IntegralCluster::PrimType &prim, \
    const SymGroup &generating_grp, \
    const jsonParser &bspecs, \
    double crystallography_tol, \
    const std::function<bool (Site)> &site_filter, \
    const typename ORBIT::SymCompareType &sym_compare, \
    INSERTER result, \
    std::ostream &status); \
  \
  template std::pair<ITERATOR, ITERATOR> orbit_branch<ITERATOR>( \
    ITERATOR begin, \
    ITERATOR end, \
    Index size); \
  \
 
#define _SPECS_IT(ORBIT) std::vector<OrbitBranchSpecs<ORBIT> >::iterator

#define _VECTOR_IT(ORBIT) std::vector<ORBIT>::iterator
#define _VECTOR_INSERTER(ORBIT) std::back_insert_iterator<std::vector<ORBIT> >

#define _SET_IT(ORBIT) std::set<ORBIT>::iterator
#define _SET_INSERTER(ORBIT) std::insert_iterator<std::set<ORBIT> >

#define CLUSTER_ORBITS_VECTOR_INST(ORBIT) CLUSTER_ORBITS_INST(_VECTOR_IT(ORBIT),_VECTOR_INSERTER(ORBIT),ORBIT,_SPECS_IT(ORBIT))
#define CLUSTER_ORBITS_SET_INST(ORBIT) CLUSTER_ORBITS_INST(_SET_IT(ORBIT),_SET_INSERTER(ORBIT),ORBIT,_SPECS_IT(ORBIT))

  CLUSTER_ORBITS_VECTOR_INST(PrimPeriodicIntegralClusterOrbit)

}
