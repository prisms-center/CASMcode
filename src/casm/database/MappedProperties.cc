#include "casm/database/MappedProperties.hh"
#include <utility>

//#include "casm/casm_io/jsonParser.hh"

namespace CASM {
  /*
  const std::string traits<ScoreMappedProperties::Method>::name = "method";

  const std::multimap<ScoreMappedProperties::Method, std::vector<std::string> > traits<ScoreMappedProperties::Method>::strval = {
    {ScoreMappedProperties::Method::deformation_cost, {"deformation_cost"} },
    {ScoreMappedProperties::Method::minimum, {"minimum"} },
    {ScoreMappedProperties::Method::maximum, {"maximum"} },
    {ScoreMappedProperties::Method::direct_selection, {"direct_selection"} }
  };

  ENUM_IO_DEF(CASM::ScoreMappedProperties::Method)
  */

  bool MappedProperties::has_scalar(std::string const &_name) const {
    auto it = global.find(_name);

    return it != global.end() && ((it->second).rows() == 1 && (it->second).cols() == 1);
  }

  double const &MappedProperties::scalar(std::string const &_name) const {
    auto it = global.find(_name);
    if(it == global.end()) {
      throw std::invalid_argument("Invalid scalar property '" + _name + "'");
    }
    assert((it->second).rows() == 1 && (it->second).cols() == 1);
    return (it->second)(0, 0);

  }

  double &MappedProperties::scalar(std::string const &_name) {
    auto it = global.find(_name);
    if(it == global.end()) {
      it = global.emplace(std::make_pair<const std::string, Eigen::MatrixXd>(_name, Eigen::MatrixXd::Zero(1, 1))).first;
    }
    assert((it->second).rows() == 1 && (it->second).cols() == 1);
    return (it->second)(0, 0);

  }
  //jsonParser &to_json(const MappedProperties &obj, jsonParser &json) {
  //json.put_obj();
  //json["from"] = obj.from;
  //json["to"] = obj.to;
  //json["unmapped"] = obj.unmapped;
  //json["mapped"] = obj.mapped;
  //return json;
  //}

  //void from_json(MappedProperties &obj, const jsonParser &json) {
  //from_json(obj.from, json["from"]);
  //from_json(obj.to, json["to"]);
  //from_json(obj.unmapped, json["unmapped"]);
  //from_json(obj.mapped, json["mapped"]);
  //}

  ScoreMappedProperties::Option::Option(Method _method, std::string _name, double _lattice_weight) :
    m_method(_method),
    m_name(std::move(_name)),
    m_lattice_weight(_lattice_weight) {
    if(_method = Method::deformation_cost) {
      if(m_lattice_weight < 0.0 || m_lattice_weight > 1.0) {
        throw std::invalid_argument("Invalid ScoreMappedProperties option, using 'deformation_cost' method, lattice_weight must be between 0.0 and 1.0");
      }
    }
    else if(m_name.empty()) {
      throw std::invalid_argument("Invalid ScoreMappedProperties option, using 'minimum', 'maximum, or 'direct_selection' method, additional string argument is required.");
    }
  }


  /// \brief Default uses minimum relaxed_energy

  ScoreMappedProperties::ScoreMappedProperties(ScoreMappedProperties::Option _opt) :
    m_opt(_opt) {

    // default to minimum "relaxed_energy"
    if(m_params.is_null() || m_params.size() == 0) {
      m_params.put_obj();
      m_params["method"] = "minimum";
      m_params["property"] = "relaxed_energy";
    }

    from_json(m_method, m_params["method"]);

    switch(m_method) {
    case Method::deformation_cost: {
      from_json(m_lattice_weight, m_params["lattice_weight"]);
      if(m_lattice_weight < 0.0 || m_lattice_weight > 1.0) {
        throw std::invalid_argument("Invalid ScoreMappedProperties lattice weight");
      }
      break;
    }
    case Method::minimum:
    case Method::maximum: {
      from_json(m_propname, m_params["property"]);
      break;
    }
    case Method::direct_selection: {
      from_json(m_direct_selection_name, m_params["name"]);
      break;
    }
    }
  }

  double ScoreMappedProperties::operator()(const MappedProperties &obj) const {
    switch(m_method) {
    case Method::minimum: {
      return obj.scalar(m_propname);
    }
    case Method::maximum: {
      return -obj.scalar(m_propname);
    }
    case Method::deformation_cost: {
      return obj.scalar("lattice_deformation_cost") * m_lattice_weight + obj.scalar("basis_deformation_cost") * (1.0 - m_lattice_weight);
    }
    case Method::direct_selection: {
      if(obj.to == m_direct_selection_name) {
        return 0.0;
      }
      else {
        return 1.0;
      }
    }
    default: {
      return std::numeric_limits<double>::max();
    }
    }
  }

  bool ScoreMappedProperties::validate(const MappedProperties &obj) const {
    switch(m_method) {
    case Method::minimum:
    case Method::maximum: {
      return obj.global.contains(m_propname);
    }
    case Method::deformation_cost: {
      return obj.global.contains("lattice_deformation_cost") &&
             obj.global.contains("basis_deformation_cost");
    }
    case Method::direct_selection: {
      return true;
    }
    default: {
      return false;
    }
    }
  }

  bool ScoreMappedProperties::operator==(const ScoreMappedProperties &B) const {
    return m_params == B.m_params;
  }

  bool ScoreMappedProperties::operator!=(const ScoreMappedProperties &B) const {
    return !(*this == B);
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
}

