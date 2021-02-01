#include "casm/clex/io/json/ConfigMapping_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/global/definitions.hh"

namespace CASM {

//--ConfigMapping::Settings------------------------------------------------

jsonParser &to_json(ConfigMapping::Settings const &_set, jsonParser &_json) {
  _json["lattice_weight"] = _set.lattice_weight;
  _json["ideal"] = _set.ideal;
  _json["strict"] = _set.strict;
  _json["primitive_only"] = _set.primitive_only;
  _json["robust"] = _set.robust;
  _json["fix_volume"] = _set.fix_volume;
  _json["fix_lattice"] = _set.fix_lattice;
  _json["k_best"] = _set.k_best;
  if (!_set.forced_lattices.empty())
    _json["forced_lattices"] = _set.forced_lattices;
  if (!_set.filter.empty()) _json["filter"] = _set.filter;
  _json["cost_tol"] = _set.cost_tol;
  _json["min_va_frac"] = _set.min_va_frac;
  _json["max_va_frac"] = _set.max_va_frac;
  _json["max_vol_change"] = _set.max_vol_change;

  return _json;
}

jsonParser const &from_json(ConfigMapping::Settings &_set,
                            jsonParser const &_json) {
  _set.set_default();

  if (_json.contains("lattice_weight"))
    _set.lattice_weight = _json["lattice_weight"].get<double>();

  if (_json.contains("ideal")) _set.ideal = _json["ideal"].get<bool>();

  if (_json.contains("strict")) _set.strict = _json["strict"].get<bool>();

  if (_json.contains("robust")) _set.robust = _json["robust"].get<bool>();

  if (_json.contains("primitive_only"))
    _set.primitive_only = _json["primitive_only"].get<bool>();

  if (_json.contains("fix_volume"))
    _set.fix_volume = _json["fix_volume"].get<bool>();

  if (_json.contains("fix_lattice"))
    _set.fix_lattice = _json["fix_lattice"].get<bool>();

  if (_json.contains("k_best")) {
    _set.k_best = _json["k_best"].get<Index>();
  }

  if (_json.contains("forced_lattices"))
    _set.forced_lattices =
        _json["forced_lattices"].get<std::vector<std::string> >();

  if (_json.contains("filter"))
    _set.filter = _json["filter"].get<std::string>();

  if (_json.contains("cost_tol"))
    _set.cost_tol = _json["cost_tol"].get<double>();

  if (_json.contains("min_va_frac"))
    _set.min_va_frac = _json["min_va_frac"].get<double>();

  if (_json.contains("max_va_frac"))
    _set.max_va_frac = _json["max_va_frac"].get<double>();

  if (_json.contains("max_vol_change"))
    _set.max_vol_change = _json["max_vol_change"].get<double>();

  return _json;
}

}  // namespace CASM
