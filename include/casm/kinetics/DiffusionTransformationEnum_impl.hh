#include "casm/kinetics/DiffusionTransformationEnum.hh"
#include "casm/kinetics/DiffusionTransformation_impl.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/OrbitGeneration.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"

namespace CASM {

  namespace Kinetics {

    /// \brief Make PrimPeriodicDiffTransOrbits from a range of PrimPeriodicOrbit<IntegralCluster>
    ///
    /// \param begin,end Range of PrimPeriodicOrbit<IntegralCluster>
    /// \param xtal_tol Tolerance for invariants comparisons
    /// \param result PrimPeriodicDiffTransOrbit output iterator
    /// \param primclex PrimClex pointer. May be nullptr if only the GenericOrbit
    ///        interface is needed, but this is probably the less likely use case
    ///        so it must be explicitly given.
    template<typename OrbitOutputIterator, typename IntegralClusterOrbitInputIterator>
    OrbitOutputIterator make_prim_periodic_diff_trans_orbits(
      IntegralClusterOrbitInputIterator begin,
      IntegralClusterOrbitInputIterator end,
      double xtal_tol,
      OrbitOutputIterator result,
      const PrimClex *primclex) {

      const Structure &prim = begin->prototype().prim();
      const auto &generating_grp = prim.factor_group();
      PrimPeriodicDiffTransSymCompare sym_compare {xtal_tol};

      OrbitGenerators<PrimPeriodicDiffTransOrbit> generators {generating_grp, sym_compare};

      for(auto it = begin; it != end; ++it) {
        DiffusionTransformationEnum e {it->prototype()};
        for(auto it = e.begin(); it != e.end(); ++it) {
          // generators.insert generates the sorted canonical form and inserts it
          generators.insert(*it);
        }
      }

      return generators.make_orbits(result, *primclex);
    }

  }

}
