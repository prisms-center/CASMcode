#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clusterography/ClusterInvariants_impl.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/symmetry/Orbit_impl.hh"

namespace CASM {

  /// \brief Generate clusters using all Site
  bool all_sites_filter(const Site &site) {
    return true;
  }

  /// \brief Generate clusters using Site with site_occupant.size() > 1
  bool alloy_sites_filter(const Site &site) {
    return site.site_occupant().size() > 1;
  }

  // --- ClusterInvariants instantiations ----------------------------------------

  // Here we are using ELEMENT as in Orbit<ELEMENT,SYMCOMPARETYPE>

#define _ELEMENT(ORBIT) typename ORBIT::Element

#define CLUSTER_INVARIANTS_INST(ELEMENT) \
  template class ClusterInvariants<ELEMENT>; \
  template bool almost_equal<ELEMENT>(const ClusterInvariants<ELEMENT> &A, const ClusterInvariants<ELEMENT> &B, double tol); \
  template bool compare<ELEMENT>(const ClusterInvariants<ELEMENT> &A, const ClusterInvariants<ELEMENT> &B, double tol); \
  template std::ostream &operator<< <ELEMENT>(std::ostream &sout, const ClusterInvariants<ELEMENT> &invariants); \

  CLUSTER_INVARIANTS_INST(IntegralCluster)


  // --- Orbit<ClusterType, SymCompareType> instantiations -----------------------

#define  CLUSTER_ORBITS_INST(ITERATOR,INSERTER,ORBIT,ELEMENT,SPECSITERTOR) \
  \
  template class ORBIT; \
  template ITERATOR find_orbit<ITERATOR,ELEMENT>(ITERATOR begin, ITERATOR end, ELEMENT e); \
  \
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

#define _ORBIT(ELEMENT,SYMCOMPARE) Orbit<ELEMENT,SYMCOMPARE>

#define CLUSTER_ORBITS_VECTOR_INST(ELEMENT,SYMCOMPARE) \
  CLUSTER_ORBITS_INST( \
    _VECTOR_IT(_ORBIT(ELEMENT,SYMCOMPARE)), \
    _VECTOR_INSERTER(_ORBIT(ELEMENT,SYMCOMPARE)), \
    _ORBIT(ELEMENT,SYMCOMPARE), \
    _ELEMENT(_ORBIT(ELEMENT,SYMCOMPARE)), \
    _SPECS_IT(_ORBIT(ELEMENT,SYMCOMPARE)))

#define CLUSTER_ORBITS_SET_INST(ELEMENT,SYMCOMPARE) \
  CLUSTER_ORBITS_INST( \
    _SET_IT(_ORBIT(ELEMENT,SYMCOMPARE)), \
    _SET_INSERTER(_ORBIT(ELEMENT,SYMCOMPARE)), \
    _ORBIT(ELEMENT,SYMCOMPARE), \
    _ELEMENT(_ORBIT(ELEMENT,SYMCOMPARE)), \
    _SPECS_IT(_ORBIT(ELEMENT,SYMCOMPARE)))

  CLUSTER_ORBITS_VECTOR_INST(IntegralCluster, LocalSymCompare<IntegralCluster>)
  CLUSTER_ORBITS_VECTOR_INST(IntegralCluster, PrimPeriodicSymCompare<IntegralCluster>)
  CLUSTER_ORBITS_VECTOR_INST(IntegralCluster, ScelPeriodicSymCompare<IntegralCluster>)

}
