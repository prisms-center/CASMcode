#ifndef UNIT_TestOrbits
#define UNIT_TestOrbits

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clusterography/ClusterOrbits.hh"

using namespace CASM;

namespace CASM {
  namespace Completer {
    class EnumOption;
  }
}

namespace test {

  struct TestPrimPeriodicIntegralClusterOrbits {

    TestPrimPeriodicIntegralClusterOrbits(
      const PrimClex &primclex,
      jsonParser _specs);

    jsonParser specs;
    std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
  };

  struct TestPrimPeriodicDiffusionTransformationOrbits : TestPrimPeriodicIntegralClusterOrbits {

    /// Make PrimPeriodicDiffusionTransformationOrbit for all
    ///   PrimPeriodicIntegralClusterOrbit branches
    TestPrimPeriodicDiffusionTransformationOrbits(
      const PrimClex &primclex,
      jsonParser _specs);

    /// Make PrimPeriodicDiffusionTransformationOrbit for a range of
    ///   PrimPeriodicIntegralClusterOrbit branches
    TestPrimPeriodicDiffusionTransformationOrbits(
      const PrimClex &primclex,
      jsonParser _specs,
      Index branch_begin,
      Index branch_end);

    /// Make PrimPeriodicDiffusionTransformationOrbit using DiffusionTransformation::run
    TestPrimPeriodicDiffusionTransformationOrbits(
      const PrimClex &primclex,
      jsonParser diff_trans_specs,
      const Completer::EnumOption &enum_opt);

    std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
  };

  struct TestLocalOrbits {
    TestLocalOrbits(
      const PrimClex &primclex,
      const Kinetics::DiffusionTransformation &diff_trans,
      jsonParser _specs);

    jsonParser specs;
    std::vector<LocalIntegralClusterOrbit> orbits;
  };
}

#endif
