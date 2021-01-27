#ifndef CASM_clex_io_json_ConfigMapping
#define CASM_clex_io_json_ConfigMapping

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/global/definitions.hh"
namespace CASM {
namespace ConfigMapping {
struct Settings;
}

jsonParser &to_json(ConfigMapping::Settings const &_set, jsonParser &_json);

jsonParser const &from_json(ConfigMapping::Settings &_set,
                            jsonParser const &_json);

}  // namespace CASM

#endif
