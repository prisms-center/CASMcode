#ifndef CASM_ClusterOrbits_json_io
#define CASM_ClusterOrbits_json_io

#include <memory>
#include <vector>

namespace CASM {

struct IntegralClusterOrbitGenerator;
class Structure;
template <typename T>
class InputParser;
class jsonParser;

/// Write custom orbit specs to JSON
jsonParser &to_json(const IntegralClusterOrbitGenerator &orbit_generator,
                    jsonParser &json);

/// Parse custom orbit specs from JSON
void parse(InputParser<std::vector<IntegralClusterOrbitGenerator>> &parser,
           Structure const &prim);

}  // namespace CASM

#endif
