#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/OrbitGeneration.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {

  namespace Kinetics {

    template<typename OrbitOutputIterator, typename IntegralClusterOrbitInputIterator>
    OrbitOutputIterator make_prim_periodic_diff_trans_orbits(
      IntegralClusterOrbitInputIterator begin,
      IntegralClusterOrbitInputIterator end,
      double xtal_tol,
      OrbitOutputIterator result) {

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

      return generators.make_orbits(result);
    }

  }

}
