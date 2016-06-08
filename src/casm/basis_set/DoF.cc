#include "casm/basis_set/DoF.hh"

// needed for to/from json... this should be fixed
#include "casm/crystallography/Molecule.hh"

namespace CASM {

  DoF::TraitsMap &DoF::_traits_map() {
    static TraitsMap _static_traits_map([](const TraitsMap::value_type & value)->std::string {
      return value.type_name();
    });
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
    m_type_name(_type_func()->type_name()),
    m_var_name(_var_name),
    m_dof_ID(_ID),
    m_ID_lock(false) {

    _traits_map()[type_name()] = *_type_func();

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
  bool DoFSet::update_IDs(const std::vector<Index> &before_IDs, const std::vector<Index> &after_IDs) {

    Index ID_ind;
    bool is_updated(false);
    for(Index i = 0; i < m_components().size(); i++) {
      // IMPORTANT: Do before_IDs.find(), NOT m_components().find() (if such a thing existed)
      ID_ind = find_index(before_IDs, m_components[i].ID());
      // Only set ID if DoF doesn't have an ID lock
      if(ID_ind < after_IDs.size() && !m_components[i].is_locked()) {
        m_components[i].set_ID(after_IDs[ID_ind]);
        // The new ID only changes the formula if the corresponding coeff is nonzero
        //if(!almost_zero(m_coeffs[i]))
        is_updated = true;
      }
    }
    return is_updated;
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

