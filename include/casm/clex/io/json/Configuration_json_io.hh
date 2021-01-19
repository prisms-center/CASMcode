#ifndef CASM_clex_Configuration_json_io
#define CASM_clex_Configuration_json_io

#include <string>

namespace CASM {

class PrimClex;
class Supercell;
class jsonParser;

template <typename T>
struct jsonConstructor;

template <>
struct jsonConstructor<Configuration> {
  /// Read Configuration from JSON
  static Configuration from_json(
      jsonParser const &json,
      std::shared_ptr<const Structure> const &shared_prim);
};

template <>
struct jsonMake<Configuration> {
  /// Read Configuration from JSON
  static std::unique_ptr<Configuration> make_from_json(
      jsonParser const &json,
      std::shared_ptr<const Structure> const &shared_prim);
};

/// Insert Configuration to JSON
jsonParser &to_json(Configuration const &configuration, jsonParser &json);

/// Read Configuration from JSON
void from_json(Configuration &configuration, jsonParser const &json);

}  // namespace CASM

#endif
