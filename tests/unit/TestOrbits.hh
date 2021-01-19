#ifndef UNIT_TestOrbits
#define UNIT_TestOrbits

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clusterography/ClusterOrbits.hh"

using namespace CASM;

namespace CASM {
namespace Completer {
class EnumOption;
}
}  // namespace CASM

namespace test {

struct TestPrimPeriodicIntegralClusterOrbits {
  TestPrimPeriodicIntegralClusterOrbits(const PrimClex &primclex,
                                        jsonParser _specs);

  jsonParser specs;
  std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
};
}  // namespace test

#endif
