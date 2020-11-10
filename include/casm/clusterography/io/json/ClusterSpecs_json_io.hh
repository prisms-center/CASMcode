#ifndef CASM_ClusterSpecs_json_io
#define CASM_ClusterSpecs_json_io

#include <memory>

namespace CASM {

  class ClusterSpecs;
  template<typename T> class InputParser;
  class LocalMaxLengthClusterSpecs;
  class PeriodicMaxLengthClusterSpecs;
  // class WithinScelMaxLengthClusterSpecs;
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

  /// \brief Parse ClusterSpecs from JSON & validate
  void parse(
    InputParser<ClusterSpecs> &parser,
    const std::shared_ptr<Structure const> &shared_prim,
    const SymGroup &super_group);

  // /// Parse ClusterSpecs from JSON for supercell orbits & validate
  // void parse(
  //   InputParser<ClusterSpecs> &parser,
  //   const std::shared_ptr<Structure const> &shared_prim,
  //   std::vector<PermuteIterator> const &super_group,
  //   SupercellSymInfo const &sym_info);

  /// \brief Write PeriodicMaxLengthClusterSpecs to JSON
  jsonParser &to_json(
    const PeriodicMaxLengthClusterSpecs &cspecs,
    jsonParser &json);

  /// \brief Write LocalMaxLengthClusterSpecs to JSON
  jsonParser &to_json(
    const LocalMaxLengthClusterSpecs &cspecs,
    jsonParser &json);

  // /// \brief Write WithinScelMaxLengthClusterSpecs to JSON
  // jsonParser &to_json(
  //   const WithinScelMaxLengthClusterSpecs &cspecs,
  //   jsonParser &json);

  /// \brief Write PeriodicMaxLengthClusterSpecs or LocalMaxLengthClusterSpecs to JSON
  jsonParser &to_json(
    const ClusterSpecs &cspecs,
    jsonParser &json);
}

#endif
