#include "casm/basis_set/DoF.hh"

// needed for to/from json... this should be fixed
#include "casm/crystallography/Molecule.hh"

namespace CASM {

  DoF::TraitsMap const &DoF::_traits_map() {
    static TraitsMap _static_traits_map([](const value_type & value)->std::string {
      return value.name();
    })
    return _static_traits_map;
  }

  //********************************************************************

  DoF::BasicTraits const &DoF::traits(std::string const &_type_name) {
    auto it = _traits_map().find(_type_name);
    if(it == _traits_map().end()) {
      throw std::runtime_error("Could not find DoF Traits for DoF type" + _type_name);
    }
    return *it;
  }

  //********************************************************************

  DoF::DoF(DoF::TypeFunc _type_func,
           std::string const &_var_name,
           Index _ID) :
    m_type_name(_type_func()->name()),
    m_var_name(_var_name),
    m_dof_ID(_ID),
    m_ID_lock(false) {

    _traits_map()[type_name()] = *_type_func();

  }

  //********************************************************************
  // ** OccupantDoF **
  /// overload for each template type to be used
  //
  //********************************************************************
  // int version
  template<> jsonParser &OccupantDoF<int>::to_json(jsonParser &json) const {
    json.put_obj();
    json["DoF_type"] = "OccupantDoF";
    json["DoF_template_type"] = "int";
    json["type_name"] = m_type_name;
    json["domain"] = m_domain;
    json["current_state"] = m_current_state;
    return json;
  }

  //********************************************************************

  jsonParser &to_json(const OccupantDoF<int> &dof, jsonParser &json) {
    return dof.to_json(json);
  }

  //********************************************************************

  void from_json(OccupantDoF<int> &dof, const jsonParser &json) {
    std::string name = json["type_name"].get<std::string>();
    std::vector<int> domain = json["domain"].get<std::vector<int> >();
    int current_state = json["current_state"].get<int>();

    dof = OccupantDoF<int>(name, domain, current_state);
  }

  //********************************************************************
  // molecule version
  template<> jsonParser &OccupantDoF<Molecule>::to_json(jsonParser &json) const {
    json.put_obj();
    json["DoF_type"] = "OccupantDoF";
    json["DoF_template_type"] = "Molecule";
    json["type_name"] = m_type_name;

    json["domain"] = m_domain;
    json["current_state"] = m_current_state;

    return json;
  }

  //********************************************************************

  jsonParser &to_json(const OccupantDoF<Molecule> &dof, jsonParser &json) {
    return dof.to_json(json);
  }

  //********************************************************************

  void from_json(OccupantDoF<Molecule> &dof, const jsonParser &json, Eigen::Matrix3d const &f2c_mat) {
    std::string name = json["type_name"].get<std::string>();

    int current_state = json["current_state"].get<int>();

    std::vector<Molecule> domain;
    from_json(domain, json["domain"], f2c_mat);

    dof = OccupantDoF<Molecule>(name, domain, current_state);

  }

  //********************************************************************
  // ** ContinuousDoF **

  jsonParser &to_json(const ContinuousDoF &dof, jsonParser &json) {
    return dof.to_json(json);
  }

  //********************************************************************

  void from_json(ContinuousDoF &dof, const jsonParser &json) {
    dof = ContinuousDoF(json["type_name"].get<std::string>(), json["min"].get<double>(), json["max"].get<double>());
  }

  //********************************************************************

  void DoFSet::from_json(jsonParser const &json) {
    /*
      if(json.contains("coordinate_space"))
      m_coordinate_space=json["coordinate_space"].get<Eigen::MatrixXd>();
    */
    json.get_if(m_excluded_occs, "excluded_occupants");

    DoF::traits(name()).from_json(*this, json);

  }

  //********************************************************************

  jsonParser &to_json(jsonParser &json) const {
    /*
    if(!m_coordinate_space.isIdentity(1e-5))
      json["coordinate_space"]=coordinate_space();
    */
    if(m_excluded_occs.count())
      json["excluded_occupants"] = m_excluded_occs;

    DoF::traits(name()).to_json(*this, json);
  }

}

