#ifndef CASM_ConfigDoF_json_io
#define CASM_ConfigDoF_json_io

#include <map>
#include <vector>

#include "casm/crystallography/DoFDecl.hh"
#include "casm/global/definitions.hh"

namespace CASM {

class ConfigDoF;
class Structure;
template <typename T>
struct jsonConstructor;
template <typename T>
struct jsonMake;
class jsonParser;

template <>
struct jsonMake<ConfigDoF> {
  static std::unique_ptr<ConfigDoF> make_from_json(const jsonParser &json,
                                                   Structure const &prim,
                                                   Index volume);
};

template <>
struct jsonConstructor<ConfigDoF> {
  static ConfigDoF from_json(jsonParser const &json, Structure const &prim,
                             Index volume);
};

jsonParser &to_json(const ConfigDoF &configdof, jsonParser &json);

void from_json(ConfigDoF &configdof, const jsonParser &json);

}  // namespace CASM

#endif
