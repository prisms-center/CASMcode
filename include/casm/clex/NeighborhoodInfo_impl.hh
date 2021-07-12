#include "casm/clex/NeighborList.hh"
#include "casm/clex/NeighborhoodInfo.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"

namespace CASM {

template <typename ClusterOrbitIterator>
NeighborhoodInfo::NeighborhoodInfo(PrimNeighborList &prim_neighbor_list,
                                   ClusterOrbitIterator begin,
                                   ClusterOrbitIterator end) {
  // construct neighborhood and "flower" function info (point corr)
  this->site_dependency_neighborhoods =
      make_site_dependency_neighborhoods(begin, end);
  this->neighbor_list_size = 0;
  this->n_point_corr = site_dependency_neighborhoods.size();
  for (auto const &neighborhood : site_dependency_neighborhoods) {
    Index key_index = prim_neighbor_list.neighbor_index(neighborhood.first);
    this->n_point_corr = max(key_index + 1, this->n_point_corr);
    for (auto const &unitcellcoord : neighborhood.second) {
      this->unitcell_neighborhood.insert(unitcellcoord.unitcell());
      this->unitcellcoord_neighborhood.insert(unitcellcoord);
      Index neighbor_index = prim_neighbor_list.neighbor_index(unitcellcoord);
      this->neighbor_list_size =
          max(neighbor_index + 1, this->neighbor_list_size);
    }
  }

  // construct vector of neighbor list sites for origin
  std::vector<xtal::UnitCellCoord> tmp_nlist_sites;
  auto it = prim_neighbor_list.begin();
  while (it != prim_neighbor_list.end()) {
    auto b_it = prim_neighbor_list.sublat_indices().begin();
    while (b_it != prim_neighbor_list.sublat_indices().end()) {
      tmp_nlist_sites.emplace_back(*b_it, *it);
      ++b_it;
    }
    ++it;
  }

  // construct point_corr_unitcellcoord
  for (Index neighbor_index = 0; neighbor_index < n_point_corr;
       ++neighbor_index) {
    this->point_corr_unitcellcoord.push_back(tmp_nlist_sites[neighbor_index]);
  }

  // find asymmetric unit indices
  for (xtal::UnitCellCoord unitcellcoord : point_corr_unitcellcoord) {
    this->point_corr_asymmetric_unit_indices.push_back(
        find_asymmetric_unit_index(unitcellcoord, begin, end));
  }
}

}  // namespace CASM
