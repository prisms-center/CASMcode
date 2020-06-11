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
  jsonParser &to_json(const IntegralCluster &clust, jsonParser &json);

  /// \brief Read from JSON
  void from_json(IntegralCluster &clust, const jsonParser &json, double xtal_tol);

  template<>
  struct jsonConstructor<IntegralCluster> {

    /// \brief Construct from JSON
    static IntegralCluster from_json(const jsonParser &json, const Structure &prim, double xtal_tol);
  };

  /// \brief Parse IntegralCluster from JSON
  void parse(
    InputParser<IntegralCluster> &parser,
    const std::shared_ptr<const Structure> &shared_prim);

  /// Write custom orbit specs to JSON
  jsonParser &to_json(const IntegralClusterOrbitGenerator &orbit_generator, jsonParser &json);

  /// Parse custom orbit specs from JSON
  void parse(
    InputParser<std::vector<IntegralClusterOrbitGenerator>> &parser,
    const std::shared_ptr<const Structure> &shared_prim);

}

#endif
