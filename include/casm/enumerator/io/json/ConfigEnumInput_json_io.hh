#ifndef CASM_enumerator_ConfigEnumInput_json_io
#define CASM_enumerator_ConfigEnumInput_json_io

#include <memory>
#include <vector>

namespace CASM {

  template<typename T> class InputParser;
  class jsonParser;
  class Structure;
  class Supercell;
  class Configuration;

  namespace DB {
    template<typename T> class Database;
  }

  /// Output ConfigEnumInput to JSON
  jsonParser &to_json(ConfigEnumInput const &config_enum_input, jsonParser &json);

  /// Read ConfigEnumInput from JSON
  template<>
  struct jsonConstructor<ConfigEnumInput> {

    static ConfigEnumInput from_json(
      const jsonParser &json,
      std::shared_ptr<Structure const> const &shared_prim,
      DB::Database<Supercell> &supercell_db);
  };

  /// Read ConfigEnumInput from JSON
  void parse(
    InputParser<ConfigEnumInput> &parser,
    std::shared_ptr<Structure const> const &shared_prim,
    DB::Database<Supercell> &supercell_db);

  /// Read std::vector<ConfigEnumInput> from JSON input, allowing queries from databases
  void from_json(
    std::vector<ConfigEnumInput> &config_enum_input,
    jsonParser const &json,
    std::shared_ptr<Structure const> shared_prim,
    DB::Database<Supercell> &supercell_db,
    DB::Database<Configuration> &configuration_db);


  /// Parse JSON to construct initial states for enumeration (as std::vector<ConfigEnumInput>)
  void parse(
    InputParser<std::vector<ConfigEnumInput>> &parser,
    std::shared_ptr<Structure const> shared_prim,
    DB::Database<Supercell> &supercell_db,
    DB::Database<Configuration> &configuration_db);

}

#endif
