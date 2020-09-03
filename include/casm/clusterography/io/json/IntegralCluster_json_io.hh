#ifndef CASM_IntegralCluster_json_io
#define CASM_IntegralCluster_json_io

#include <memory>

namespace CASM {

  class IntegralCluster;
  struct IntegralClusterOrbitGenerator;
  class Structure;
  template<typename T> class InputParser;
  class jsonParser;

  /// \brief Write IntegralCluster to JSON object
  jsonParser &to_json(IntegralCluster const &clust, jsonParser &json);

  /// \brief Read from JSON
  void from_json(IntegralCluster &clust, jsonParser const &json);

  template<>
  struct jsonConstructor<IntegralCluster> {

    /// \brief Construct from JSON
    static IntegralCluster from_json(jsonParser const &json, Structure const &prim);
  };

  /// \brief Parse IntegralCluster from JSON
  void parse(InputParser<IntegralCluster> &parser, Structure const &prim);

}

#endif
