#include "casm/clex/MappedProperties.hh"

#include <utility>

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"

namespace CASM {
/*
const std::string traits<ScoreMappedProperties::Method>::name = "method";


const std::multimap<ScoreMappedProperties::Method, std::vector<std::string> >
traits<ScoreMappedProperties::Method>::strval = {
  {ScoreMappedProperties::Method::deformation_cost, {"deformation_cost"} },
  {ScoreMappedProperties::Method::minimum, {"minimum"} },
  {ScoreMappedProperties::Method::maximum, {"maximum"} },
  {ScoreMappedProperties::Method::direct_selection, {"direct_selection"} }
};
*/
/*
ENUM_IO_DEF(CASM::ScoreMappedProperties::Method)
*/

bool MappedProperties::has_scalar(std::string const &_name) const {
  auto it = global.find(_name);

  return it != global.end() &&
         ((it->second).rows() == 1 && (it->second).cols() == 1);
}

double const &MappedProperties::scalar(std::string const &_name) const {
  auto it = global.find(_name);
  if (it == global.end()) {
    throw std::invalid_argument("Invalid scalar property '" + _name + "'");
  }
  assert((it->second).rows() == 1 && (it->second).cols() == 1);
  return (it->second)(0, 0);
}

double &MappedProperties::scalar(std::string const &_name) {
  auto it = global.find(_name);
  if (it == global.end()) {
    it = global
             .emplace(std::make_pair(_name, Eigen::MatrixXd::Zero(1, 1).eval()))
             .first;
  }
  assert((it->second).rows() == 1 && (it->second).cols() == 1);
  return (it->second)(0, 0);
}
jsonParser &to_json(const MappedProperties &obj, jsonParser &json) {
  json.put_obj();
  json["origin"] = obj.origin;
  json["init"] = obj.init_config;
  json["to"] = obj.to;
  if (!obj.file_data.empty()) {
    json["path"] = obj.file_data.path();
    json["timestamp"] = obj.file_data.timestamp();
  }
  json["global"] = obj.global;
  json["site"] = obj.site;

  return json;
}

jsonParser const &from_json(MappedProperties &obj, const jsonParser &json) {
  from_json(obj.origin, json["origin"]);
  from_json(obj.init_config, json["init"]);
  from_json(obj.to, json["to"]);
  from_json(obj.global, json["global"]);
  from_json(obj.site, json["site"]);
  if (json.contains("path"))
    obj.file_data = FileData(json["path"].get<std::string>(),
                             json["timestamp"].get<std::time_t>());
  return json;
}

ScoreMappedProperties::Option::Option(Method _m, std::string _name,
                                      double _lattice_weight)
    : method(_m), name(std::move(_name)), lattice_weight(_lattice_weight) {
  if (_m == Method::minimum || _m == Method::maximum) {
    if (name.empty())
      throw std::invalid_argument(
          "Scoring method '" + method_name(_m) +
          "' requires additional string field \"property\".");
  } else if (_m == Method::direct_selection) {
    if (name.empty())
      throw std::invalid_argument(
          "Scoring method '" + method_name(_m) +
          "' requires additional string field \"name\".");
  } else if (_m == Method::deformation_cost) {
    if (lattice_weight < 0. || lattice_weight > 1.)
      throw std::invalid_argument("Scoring method '" + method_name(_m) +
                                  "' requires numerical parameter "
                                  "\"lattice_weight\" between 0.0 and 1.0.");
  }
}

/// \brief Default uses minimum energy

ScoreMappedProperties::ScoreMappedProperties(ScoreMappedProperties::Option _opt)
    : m_opt(_opt) {}

double ScoreMappedProperties::operator()(const MappedProperties &obj) const {
  switch (m_opt.method) {
    case Method::minimum: {
      if (!obj.has_scalar(m_opt.name)) {
        return std::numeric_limits<double>::max();
      }
      return obj.scalar(m_opt.name);
    }
    case Method::maximum: {
      if (!obj.has_scalar(m_opt.name)) {
        return std::numeric_limits<double>::max();
      }
      return -obj.scalar(m_opt.name);
    }
    case Method::deformation_cost: {
      if (!obj.has_scalar("lattice_deformation_cost") ||
          !obj.has_scalar("atomic_deformation_cost")) {
        return std::numeric_limits<double>::max();
      }
      return obj.scalar("lattice_deformation_cost") * m_opt.lattice_weight +
             obj.scalar("atomic_deformation_cost") *
                 (1.0 - m_opt.lattice_weight);
    }
    case Method::direct_selection: {
      if (obj.to == m_opt.name) {
        return 0.0;
      } else {
        return 1.0;
      }
    }
    default: {
      return std::numeric_limits<double>::max();
    }
  }
}

bool ScoreMappedProperties::validate(const MappedProperties &obj) const {
  switch (m_opt.method) {
    case Method::minimum:
    case Method::maximum:
      return obj.global.count(m_opt.name);
    case Method::deformation_cost:
      return obj.has_scalar("lattice_deformation_cost") &&
             obj.has_scalar("atomic_deformation_cost");
    case Method::direct_selection:
      return true;
    default:
      return false;
  }
}

bool ScoreMappedProperties::operator==(const ScoreMappedProperties &B) const {
  return m_opt.method == B.m_opt.method && m_opt.name == B.m_opt.name;
}

bool ScoreMappedProperties::operator!=(const ScoreMappedProperties &B) const {
  return !(*this == B);
}

jsonParser &to_json(const ScoreMappedProperties &score, jsonParser &json) {
  switch (score.option().method) {
    case ScoreMappedProperties::Method::minimum:
      json["method"] = "minimum";
      json["property"] = score.option().name;
      break;
    case ScoreMappedProperties::Method::maximum:
      json["method"] = "maximum";
      json["property"] = score.option().name;
      break;
    case ScoreMappedProperties::Method::deformation_cost:
      json["method"] = "deformation_cost";
      break;
    case ScoreMappedProperties::Method::direct_selection:
      json["method"] = "direct_selection";
      json["name"] = score.option().name;
      break;
  }

  return json;
}

jsonParser const &from_json(ScoreMappedProperties &score,
                            const jsonParser &json) {
  ScoreMappedProperties::Option opt;
  if (json.contains("method")) {
    std::string method = json["method"].get<std::string>();
    if (method == "minimum") {
      opt.method = ScoreMappedProperties::Method::minimum;
      opt.name = json["property"].get<std::string>();
    } else if (method == "maximum") {
      opt.method = ScoreMappedProperties::Method::maximum;
      opt.name = json["property"].get<std::string>();
    } else if (method == "deformation_cost") {
      opt.method = ScoreMappedProperties::Method::deformation_cost;
      opt.name = "";
    } else if (method == "direct_selection") {
      opt.method = ScoreMappedProperties::Method::direct_selection;
      opt.name = json["name"].get<std::string>();
    }
  }
  score = ScoreMappedProperties(opt);
  return json;
}

/*
  const jsonParser &ScoreMappedProperties::params() const {
  return m_params;
  }*/

/*
  jsonParser &to_json(const ScoreMappedProperties &score, jsonParser &json) {
  json = score.params();
  return json;
  }

  void from_json(ScoreMappedProperties &score, const jsonParser &json) {
  score = ScoreMappedProperties(json);
  }
*/
}  // namespace CASM
