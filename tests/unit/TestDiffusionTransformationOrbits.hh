#ifndef UNIT_TestDiffusionTransformationOrbits
#define UNIT_TestDiffusionTransformationOrbits

using namespace CASM;

namespace test {

  struct TestDiffusionTransformationOrbits {

    TestDiffusionTransformationOrbits(
      const PrimClex &primclex,
      const jsonParser &bspecs,
      Log &log = null_log()) {

      // Make PrimPeriodicIntegralClusterOrbit
      make_prim_periodic_orbits(
        primclex.prim(),
        bspecs,
        alloy_sites_filter,
        primclex.crystallography_tol(),
        std::back_inserter(orbits),
        log);
    }

    std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
  };
}

#endif
