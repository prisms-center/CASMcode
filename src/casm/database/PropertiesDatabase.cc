#include "casm/database/PropertiesDatabase.hh"

namespace CASM {

  const std::string traits<DB::ScoreMappedProperties::Method>::name = "method";

  const std::multimap<DB::ScoreMappedProperties::Method, std::vector<std::string> > traits<DB::ScoreMappedProperties::Method>::strval = {
    {DB::ScoreMappedProperties::Method::deformation_cost, {"deformation_cost"} },
    {DB::ScoreMappedProperties::Method::minimum, {"minimum"} },
    {DB::ScoreMappedProperties::Method::maximum, {"maximum"} },
    {DB::ScoreMappedProperties::Method::direct_selection, {"direct_selection"} }
  };

}

namespace CASM {
  namespace DB {

    jsonParser &to_json(const MappedProperties &obj, jsonParser &json) {
      json.put_obj();
      json["from"] = obj.from;
      json["to"] = obj.to;
      json["unmapped"] = obj.unmapped;
      json["mapped"] = obj.mapped;
      return json;
    }

    void from_json(MappedProperties &obj, const jsonParser &json) {
      from_json(obj.from, json["from"]);
      from_json(obj.to, json["to"]);
      from_json(obj.unmapped, json["unmapped"]);
      from_json(obj.mapped, json["mapped"]);
    }

    /// \brief Default uses minimum relaxed_energy
    ScoreMappedProperties::ScoreMappedProperties() :
      ScoreMappedProperties(jsonParser()) {}

    ScoreMappedProperties::ScoreMappedProperties(const jsonParser &_params) :
      m_params(_params) {

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
        return obj.mapped[m_propname].get<double>();
      }
      case Method::maximum: {
        return -1.0 * obj.mapped[m_propname].get<double>();
      }
      case Method::deformation_cost: {
        double ld = obj.mapped["lattice_deformation_cost"].get<double>();
        double bd = obj.mapped["basis_deformation_cost"].get<double>();
        return ld * m_lattice_weight + bd * (1.0 - m_lattice_weight);
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
        return obj.mapped.contains(m_propname);
      }
      case Method::deformation_cost: {
        return obj.mapped.contains("lattice_deformation_cost") &&
               obj.mapped.contains("basis_deformation_cost");
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

    const jsonParser &ScoreMappedProperties::params() const {
      return m_params;
    }

    jsonParser &to_json(const ScoreMappedProperties &score, jsonParser &json) {
      json = score.params();
      return json;
    }

    void from_json(ScoreMappedProperties &score, const jsonParser &json) {
      score = ScoreMappedProperties(json);
    }

    /// \brief Compare mapped properties 'from_A' and 'from_B', preferring self-mapped results
    bool PropertiesDatabase::Compare::operator()(
      const std::string &from_A,
      const std::string &from_B) const {

      if(from_A == from_B) {
        return false;
      }
      if(from_A == m_to) {
        return true;
      }
      if(from_B == m_to) {
        return false;
      }
      return m_score(*m_map->find_via_from(from_A)) < m_score(*m_map->find_via_from(from_B));
    }

    /// \brief Insert data
    std::pair<PropertiesDatabase::iterator, bool> PropertiesDatabase::insert(
      const MappedProperties &value) {

      // insert data
      auto res = _insert(value);
      if(!res.second) {
        return res;
      }

      // insert 'to' -> 'from' link
      auto tset = relaxed_from_all(value.to);
      tset.insert(value.from);
      _set_relaxed_from_all(value.to, tset);

      return res;
    }

    /// \brief Erase the data 'from' from_configname
    PropertiesDatabase::iterator PropertiesDatabase::erase(iterator pos) {

      auto tset = relaxed_from_all(pos->to);
      tset.erase(pos->from);
      _set_relaxed_from_all(pos->to, tset);

      return _erase(pos);
    }

  }
}

