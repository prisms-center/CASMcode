#ifndef CASM_enumerator_ConfigEnumInput_json_io
#define CASM_enumerator_ConfigEnumInput_json_io

#include <memory>
#include <vector>

namespace CASM {

namespace xtal {
class ScelEnumProps;
}

template <typename T>
class InputParser;
class jsonParser;
template <typename T>
struct jsonConstructor;
class Configuration;
class ConfigEnumInput;
class PrimClex;
class Structure;
class Supercell;

namespace DB {
template <typename T>
class Database;
}

/// Output ConfigEnumInput to JSON
jsonParser &to_json(ConfigEnumInput const &config_enum_input, jsonParser &json);

/// Read ConfigEnumInput from JSON
template <>
struct jsonConstructor<ConfigEnumInput> {
  static ConfigEnumInput from_json(
      const jsonParser &json,
      std::shared_ptr<Structure const> const &shared_prim);
};

/// Read ConfigEnumInput from JSON
void parse(InputParser<ConfigEnumInput> &parser,
           std::shared_ptr<Structure const> const &shared_prim);

/// Read std::map<std::string, ConfigEnumInput> from JSON input, allowing
/// queries from databases
void from_json(
    std::vector<std::pair<std::string, ConfigEnumInput>> &config_enum_input,
    jsonParser const &json, std::shared_ptr<Structure const> shared_prim,
    PrimClex const *primclex, DB::Database<Supercell> &supercell_db,
    DB::Database<Configuration> &configuration_db);

/// Make a ScelEnumProps object from JSON input and support
/// "unit_cell"==supercell name
void parse(InputParser<xtal::ScelEnumProps> &parser,
           DB::Database<Supercell> &supercell_db);

/// A string describing the JSON format for parsing named ConfigEnumInput
std::string parse_ConfigEnumInput_desc();

/// Parse JSON to construct initial states for enumeration (as
/// std::map<std::string, ConfigEnumInput>)
void parse(
    InputParser<std::vector<std::pair<std::string, ConfigEnumInput>>> &parser,
    std::shared_ptr<Structure const> shared_prim, PrimClex const *primclex,
    DB::Database<Supercell> &supercell_db,
    DB::Database<Configuration> &configuration_db);

}  // namespace CASM

#endif
