#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/symmetry/PermuteIterator.hh"

namespace CASM {

  IntegralClusterOrbitGenerator::IntegralClusterOrbitGenerator(
    IntegralCluster const &_prototype,
    bool _include_subclusters):
    prototype(_prototype), include_subclusters(_include_subclusters) {}

  /// Construct "within_scel" orbits from "prim_periodic" orbits
  ///
  /// The "within_scel" orbits take into account aliasing (when clusters are equivalent due to
  /// periodic boundary conditions mapping sites back "within" the supercell) when checking for
  /// cluster equivalence.
  ///
  /// This function takes clusters in the `prim_periodic_orbits`, copies and translates them
  /// throughout the supercell defined by `shared_prim` and `transformation_matrix_to_super`, and
  /// then keeps only the unique clusters under application of the `generating_group`. These
  /// generating elements are used to construct the "within_scel" orbits.
  ///
  std::vector<WithinScelIntegralClusterOrbit> make_within_scel_orbits_from_prim_periodic(
    std::shared_ptr<Structure const> const &shared_prim,
    Eigen::Matrix3l const &transformation_matrix_to_super,
    std::vector<PermuteIterator> const &generating_group,
    std::vector<PrimPeriodicIntegralClusterOrbit> const &prim_periodic_orbits) {

    // Create `generators` constructor arguments
    auto const &T = transformation_matrix_to_super;
    xtal::Lattice super_lattice = make_superlattice(shared_prim->lattice(), T);
    SymGroup _generating_group = make_sym_group(generating_group, super_lattice);
    double tol = shared_prim->lattice().tol();
    WithinScelSymCompare<IntegralCluster> sym_compare {shared_prim, T, tol};

    // Inserting clusters into `generators` keeps only the symmetrically
    // unique clusters as defined by "within_scel" sym_compare
    OrbitGenerators<WithinScelIntegralClusterOrbit> generators {_generating_group, sym_compare};

    // Try inserting all clusters in the supercell
    auto unitcells = xtal::make_lattice_points(T);
    for(const auto &orbit : prim_periodic_orbits) {
      for(const auto &equiv : orbit) {
        for(const auto &unitcell : unitcells) {
          generators.insert(equiv + unitcell);
        }
      }
    }

    // The clusters that exist in `generators` are the symmetrically
    // unique generating elements for the "within scel" orbits
    std::vector<WithinScelIntegralClusterOrbit> result;
    generators.make_orbits(std::back_inserter(result));
    return result;
  }
}
