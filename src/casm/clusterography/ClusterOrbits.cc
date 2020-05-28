#include "casm/clusterography/ClusterOrbits_impl.hh"

namespace CASM {

  IntegralClusterOrbitGenerator::IntegralClusterOrbitGenerator(
    IntegralCluster const &_prototype,
    bool _include_subclusters):
    prototype(_prototype), include_subclusters(_include_subclusters) {}

}

// namespace CASM {
//   // --- Orbit<ClusterType, SymCompareType> instantiations -----------------------
//
// #define  CLUSTER_ORBITS_INST(ITERATOR,INSERTER,ORBIT,ELEMENT,SPECSITERTOR) \
//   \
//   template class ORBIT; \
//   template ITERATOR find_orbit<ITERATOR,ELEMENT>(ITERATOR begin, ITERATOR end, ELEMENT e); \
//   \
//   template INSERTER make_asymmetric_unit<ORBIT,INSERTER>( \
//     const OrbitBranchSpecs<ORBIT> &specs, \
//     INSERTER result, \
//     std::ostream &status); \
//   \
//   template INSERTER make_next_orbitbranch<ORBIT,ITERATOR,INSERTER>( \
//     ITERATOR begin, \
//     ITERATOR end, \
//     const OrbitBranchSpecs<ORBIT> &specs, \
//     INSERTER result, \
//     std::ostream &status); \
//   \
//   template INSERTER make_orbits<SPECSITERTOR,INSERTER>( \
//     SPECSITERTOR begin, \
//     SPECSITERTOR end, \
//     const std::vector<IntegralCluster> &custom_generators, \
//     INSERTER result, \
//     std::ostream &status); \
//   \
//   template std::pair<ITERATOR, ITERATOR> orbit_branch<ITERATOR>( \
//     ITERATOR begin, \
//     ITERATOR end, \
//     Index size); \
//   \
//
// #define  PRIM_PERIODIC_CLUSTER_ORBITS_INST(INSERTER) \
//   \
//   template INSERTER make_prim_periodic_asymmetric_unit<INSERTER>( \
//     const IntegralCluster::PrimType_ptr prim_ptr, \
//     const std::function<bool (Site)> &site_filter, \
//     double xtal_tol, \
//     INSERTER result, \
//     std::ostream &status); \
//   \
//   template INSERTER make_prim_periodic_orbits<INSERTER>( \
//     const IntegralCluster::PrimType_ptr prim_ptr, \
//     const std::vector<double> &max_length, \
//     const std::vector<IntegralCluster> &custom_generators, \
//     const std::function<bool (Site)> &site_filter, \
//     double xtal_tol, \
//     INSERTER result, \
//     std::ostream &status); \
//   \
//   template INSERTER make_prim_periodic_orbits<INSERTER>( \
//     const IntegralCluster::PrimType_ptr prim_ptr, \
//     const jsonParser &bspecs, \
//     const std::function<bool (Site)> &site_filter, \
//     double xtal_tol, \
//     INSERTER result, \
//     std::ostream &status); \
//   \
//
// #define LOCAL_CLUSTER_ORBITS_INST(INSERTER,SYMCOMPARE) \
//   \
//   template INSERTER make_local_orbits<INSERTER,SYMCOMPARE>( \
//     const Kinetics::DiffusionTransformation &diff_trans, \
//     const SymGroup &generating_group, \
//     const SYMCOMPARE &sym_compare, \
//     const std::vector<double> &cutoff_radius, \
//     const std::vector<double> &max_length, \
//     const std::vector<IntegralCluster> &custom_generators, \
//     const std::function<bool (Site)> &site_filter, \
//     double xtal_tol, \
//     INSERTER result, \
//     std::ostream &status); \
//   \
//   template INSERTER make_local_orbits<INSERTER,SYMCOMPARE>( \
//     const Kinetics::DiffusionTransformation &diff_trans, \
//     const SymGroup &generating_group, \
//     const SYMCOMPARE &sym_compare, \
//     const jsonParser &bspecs, \
//     const std::function<bool (Site)> &site_filter, \
//     double xtal_tol, \
//     INSERTER result, \
//     std::ostream &status); \
//   \
//
// #define _SPECS_IT(ORBIT) std::vector<OrbitBranchSpecs<ORBIT> >::iterator
//
// #define _VECTOR_IT(ORBIT) std::vector<ORBIT>::iterator
// #define _VECTOR_INSERTER(ORBIT) std::back_insert_iterator<std::vector<ORBIT> >
//
// #define _ORBIT(SYMCOMPARE) GenericOrbit<SYMCOMPARE>
//
// #define CLUSTER_ORBITS_VECTOR_INST(SYMCOMPARE) \
//   CLUSTER_ORBITS_INST( \
//     _VECTOR_IT(_ORBIT(SYMCOMPARE)), \
//     _VECTOR_INSERTER(_ORBIT(SYMCOMPARE)), \
//     _ORBIT(SYMCOMPARE), \
//     _ELEMENT(_ORBIT(SYMCOMPARE)), \
//     _SPECS_IT(_ORBIT(SYMCOMPARE)))
//
//
//   CLUSTER_ORBITS_VECTOR_INST(LocalSymCompare<IntegralCluster>)
//   CLUSTER_ORBITS_VECTOR_INST(PrimPeriodicSymCompare<IntegralCluster>)
//   CLUSTER_ORBITS_VECTOR_INST(ScelPeriodicSymCompare<IntegralCluster>)
//   CLUSTER_ORBITS_VECTOR_INST(WithinScelSymCompare<IntegralCluster>)
//
// #define PRIM_PERIODIC_CLUSTER_ORBITS_VECTOR_INST(SYMCOMPARE) \
//   PRIM_PERIODIC_CLUSTER_ORBITS_INST( \
//     _VECTOR_INSERTER(_ORBIT(SYMCOMPARE)))
//
//   PRIM_PERIODIC_CLUSTER_ORBITS_VECTOR_INST(PrimPeriodicSymCompare<IntegralCluster>)
//
// #define LOCAL_CLUSTER_ORBITS_VECTOR_INST(SYMCOMPARE) \
//   LOCAL_CLUSTER_ORBITS_INST( \
//     _VECTOR_INSERTER(_ORBIT(SYMCOMPARE)), \
//     SYMCOMPARE)
//
//   LOCAL_CLUSTER_ORBITS_VECTOR_INST(LocalSymCompare<IntegralCluster>)
//   LOCAL_CLUSTER_ORBITS_VECTOR_INST(ScelPeriodicSymCompare<IntegralCluster>)
//   LOCAL_CLUSTER_ORBITS_VECTOR_INST(WithinScelSymCompare<IntegralCluster>)
//
// }
