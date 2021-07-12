#include "casm/casm_io/Log.hh"
#include "casm/clex/NeighborhoodInfo_impl.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"

namespace CASM {

namespace {
struct MakeNeighborhoodInfo {
  MakeNeighborhoodInfo(PrimNeighborList &_prim_neighbor_list,
                       std::unique_ptr<NeighborhoodInfo> &_neighborhood_info)
      : prim_neighbor_list(_prim_neighbor_list),
        neighborhood_info(_neighborhood_info) {}

  PrimNeighborList &prim_neighbor_list;
  std::unique_ptr<NeighborhoodInfo> &neighborhood_info;

  template <typename OrbitVecType>
  void operator()(OrbitVecType const &orbits) {
    neighborhood_info = notstd::make_unique<NeighborhoodInfo>(
        prim_neighbor_list, orbits.begin(), orbits.end());
  }
};
}  // namespace

std::unique_ptr<NeighborhoodInfo> make_neighborhood_info(
    ClusterSpecs const &cluster_specs, PrimNeighborList &prim_neighbor_list) {
  std::unique_ptr<NeighborhoodInfo> neighborhood_info;
  MakeNeighborhoodInfo f{prim_neighbor_list, neighborhood_info};
  for_all_orbits(cluster_specs, CASM::log(), f);
  return neighborhood_info;
}

std::unique_ptr<NeighborhoodInfo> make_neighborhood_info(
    ClusterSpecs const &cluster_specs,
    std::vector<IntegralCluster> const &generating_elements,
    PrimNeighborList &prim_neighbor_list) {
  std::unique_ptr<NeighborhoodInfo> neighborhood_info;
  MakeNeighborhoodInfo f{prim_neighbor_list, neighborhood_info};
  for_all_orbits(cluster_specs, generating_elements, f);
  return neighborhood_info;
}

}  // namespace CASM
