#include "casm/clex/ClexBasisInfo.hh"
#include <memory>

#include "casm/casm_io/Log.hh"
#include "casm/clex/ClexBasis_impl.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/clusterography/io/OrbitPrinter_impl.hh"

namespace CASM
{

namespace
{

/// Implementation for `make_clex_basis_info`
class MakeClexBasisInfo
{
public:
    /// Required input & resulting `clex_basis_info` to be populated
    MakeClexBasisInfo(std::shared_ptr<Structure const> _shared_prim, ClexBasisSpecs const& _basis_set_specs,
                      PrimNeighborList& _prim_neighbor_list, ClexBasisInfo& _clex_basis_info)
        : shared_prim(_shared_prim),
          basis_set_specs(_basis_set_specs),
          prim_neighbor_list(_prim_neighbor_list),
          clex_basis_info(_clex_basis_info)
    {
    }

    /// Gets called by `for_all_orbits`
    template <typename OrbitVecType> void operator()(OrbitVecType const& orbits) const
    {
        _make_neighborhood_info(orbits);
        _make_basis_function_info(orbits);
    }

    std::shared_ptr<Structure const> shared_prim;
    ClexBasisSpecs basis_set_specs;
    PrimNeighborList& prim_neighbor_list;
    ClexBasisInfo& clex_basis_info;

private:
    template <typename OrbitVecType> void _make_neighborhood_info(OrbitVecType const& orbits) const
    {
        // construct neighborhood and "flower" function info (point corr)
        auto site_dependency_neighborhoods = make_site_dependency_neighborhoods(orbits.begin(), orbits.end());
        Index neighborhood_size = 0;
        Index n_point_corr = site_dependency_neighborhoods.size();
        for (auto const& neighborhood : site_dependency_neighborhoods)
        {
            Index key_index = prim_neighbor_list.neighbor_index(neighborhood.first);
            n_point_corr = max(key_index + 1, n_point_corr);
            for (auto const& unitcellcoord : neighborhood.second)
            {
                Index neighbor_index = prim_neighbor_list.neighbor_index(unitcellcoord);
                neighborhood_size = max(neighbor_index + 1, neighborhood_size);
            }
        }
        clex_basis_info.site_dependency_neighborhoods = site_dependency_neighborhoods;
        clex_basis_info.neighborhood_size = neighborhood_size;
        clex_basis_info.n_point_corr = n_point_corr;
    }

    template <typename OrbitVecType> void _make_basis_function_info(OrbitVecType const& orbits) const
    {
        // construct ClexBasis
        ParsingDictionary<DoFType::Traits> const* dof_dict = &DoFType::traits_dict();
        ClexBasis clex_basis{shared_prim, basis_set_specs, dof_dict};
        clex_basis.generate(orbits.begin(), orbits.end());

        // collect basis function info
        Index orbit_index = 0;
        Index function_index = 0;
        for (auto const& orbit : orbits)
        {
            std::vector<Index> invariant_group_indices;
            for (SymOp const& op : orbit.equivalence_map()[0])
            {
                invariant_group_indices.push_back(op.index());
            }

            BasisSet const& bset_prototype = clex_basis.bset_orbit(orbit_index)[0];
            for (Index i = 0; i < bset_prototype.size(); ++i)
            {
                clex_basis_info.basis_function_info.emplace_back(orbit_index, function_index, orbit.prototype(),
                                                                 orbit.size(), invariant_group_indices);
                function_index++;
            }
            orbit_index++;
        }
    }
};

} // namespace

/// Make ClexBasisInfo
ClexBasisInfo make_clex_basis_info(std::shared_ptr<Structure const> const& shared_prim,
                                   ClexBasisSpecs const& basis_set_specs, PrimNeighborList& prim_neighbor_list)
{
    ClexBasisInfo clex_basis_info;
    MakeClexBasisInfo f{shared_prim, basis_set_specs, prim_neighbor_list, clex_basis_info};
    for_all_orbits(*basis_set_specs.cluster_specs, log(), f);
    return clex_basis_info;
}
} // namespace CASM
