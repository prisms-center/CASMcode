#ifndef CASM_ClusterSpecs_json_io
#define CASM_ClusterSpecs_json_io

#include <memory>

namespace CASM {

  class ClusterSpecs;
  template<typename T> class InputParser;
  class LocalMaxLengthClusterSpecs;
  class PeriodicMaxLengthClusterSpecs;
  class PrimClex;
  class SymGroup;
  class jsonParser;
  class Structure;


  /// Parse PeriodicMaxLengthClusterSpecs from JSON
  void parse(
    InputParser<PeriodicMaxLengthClusterSpecs> &parser,
    const std::shared_ptr<Structure const> &shared_prim,
    const SymGroup &super_group);

  /// Parse LocalMaxLengthClusterSpecs from JSON
  void parse(
    InputParser<LocalMaxLengthClusterSpecs> &parser,
    const std::shared_ptr<Structure const> &shared_prim,
    const SymGroup &super_group);

  /// \brief Parse PeriodicMaxLengthClusterSpecs or LocalMaxLengthClusterSpecs from JSON & validate
  void parse(
    InputParser<ClusterSpecs> &parser,
    const std::shared_ptr<Structure const> &shared_prim);

  /// \brief Parse PeriodicMaxLengthClusterSpecs or LocalMaxLengthClusterSpecs from JSON &
  ///        validate with specified generating_group
  void parse(
    InputParser<ClusterSpecs> &parser,
    const std::shared_ptr<Structure const> &shared_prim,
    const SymGroup &super_group);

  /// \brief Write PeriodicMaxLengthClusterSpecs to JSON
  jsonParser &to_json(
    const PeriodicMaxLengthClusterSpecs &cspecs,
    jsonParser &json);

  /// \brief Write LocalMaxLengthClusterSpecs to JSON
  jsonParser &to_json(
    const LocalMaxLengthClusterSpecs &cspecs,
    jsonParser &json);

  /// \brief Write PeriodicMaxLengthClusterSpecs or LocalMaxLengthClusterSpecs to JSON
  jsonParser &to_json(
    const ClusterSpecs &cspecs,
    jsonParser &json);
}

#endif
