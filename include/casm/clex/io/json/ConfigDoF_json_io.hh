#ifndef CASM_ConfigDoF_json_io
#define CASM_ConfigDoF_json_io

#include <map>
#include <vector>
#include "casm/global/definitions.hh"
#include "casm/crystallography/DoFDecl.hh"

namespace CASM {

  class ConfigDoF;
  struct DoFSetInfo;
  class LocalDiscreteConfigDoFValues;
  class LocalContinuousConfigDoFValues;
  class GlobalContinuousConfigDoFValues;
  class Structure;
  class SymGroupRepID;
  template<typename T> struct jsonMake;
  template<typename T> struct jsonConstructor;
  class jsonParser;

  jsonParser &to_json(LocalDiscreteConfigDoFValues const &_values, jsonParser &_json);
  void from_json(LocalDiscreteConfigDoFValues &_values, jsonParser const &_json);

  jsonParser &to_json(LocalContinuousConfigDoFValues const &_values, jsonParser &_json);
  void from_json(LocalContinuousConfigDoFValues &_values, jsonParser const &_json);

  jsonParser &to_json(GlobalContinuousConfigDoFValues const &_values, jsonParser &_json);
  void from_json(GlobalContinuousConfigDoFValues &_values, jsonParser const &_json);


  template<>
  struct jsonMake<ConfigDoF> {

    static std::unique_ptr<ConfigDoF> make_from_json(const jsonParser &json,
                                                     Structure const &prim);
  };

  template<>
  struct jsonConstructor<ConfigDoF> {

    static ConfigDoF from_json(jsonParser const &json, Structure const &prim);

  };

  jsonParser &to_json(const ConfigDoF &configdof, jsonParser &json);

  void from_json(ConfigDoF &configdof, const jsonParser &json);

}

#endif
