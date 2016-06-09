#include "casm/basis_set/DoF.hh"

// needed for to/from json... this should be fixed
#include "casm/crystallography/Molecule.hh"

namespace CASM {
  namespace DoF_impl {
    /// \brief implements json parsing of a specialized DoF.
    // In future, we may need to add another inheritance layer to handle DiscreteDoF types
    void BasicTraits::from_json(ContinuousDoF &_dof, jsonParser const &json) const {
      return;
    }

    /// \brief implements json parsing of a specialized DoF.
    // In future, we may need to add another inheritance layer to handle DiscreteDoF types
    void BasicTraits::to_json(ContinuousDoF const &_dof, jsonParser &json) const {
      return;
    }

    /// \brief implements json parsing of a specialized DoFSet.
    void BasicTraits::from_json(DoFSet &_dof, jsonParser const &json) const {
      return;
    }

    /// \brief implements json parsing of a specialized DoFSet.
    void BasicTraits::to_json(DoFSet const &_dof, jsonParser &json) const {
      return;
    }
  }

  DoF::TraitsMap &DoF::_traits_map() {
    static TraitsMap _static_traits_map([](const TraitsMap::value_type & value)->std::string {
      return value.type_name();
    },
    DoF_impl::traits2cloneable_ptr);
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

  DoF::DoF(DoF::BasicTraits const &_traits,
           std::string const &_var_name,
           Index _ID) :
    m_type_name(_traits.type_name()),
    m_var_name(_var_name),
    m_dof_ID(_ID),
    m_ID_lock(false) {

    _traits_map().insert(_traits);

  }

  //********************************************************************
  // ** ContinuousDoF **

  jsonParser &to_json(const ContinuousDoF &dof, jsonParser &json) {
    return dof.to_json(json);
  }

  //********************************************************************

  void from_json(ContinuousDoF &dof, const jsonParser &json) {
    DoF::traits(dof.type_name()).from_json(dof, json);
  }

  //********************************************************************

  ContinuousDoF jsonConstructor<ContinuousDoF>::from_json(const jsonParser &json, DoF::BasicTraits const &_traits) {
    ContinuousDoF tdof(_traits);
    CASM::from_json(tdof, json);
    return tdof;
  }

  //********************************************************************

  bool DoFSet::identical(DoFSet const &rhs)const {
    if(!std::equal(begin(), end(), rhs.begin(), compare_no_value)) {
      return false;
    }

    return almost_equal(m_coordinate_space, rhs.m_coordinate_space);
  }

  //********************************************************************

  bool DoFSet::update_IDs(const std::vector<Index> &before_IDs, const std::vector<Index> &after_IDs) {

    Index ID_ind;
    bool is_updated(false);
    for(Index i = 0; i < m_components.size(); i++) {
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

    DoF::traits(type_name()).from_json(*this, json);

  }

  //********************************************************************

  jsonParser &DoFSet::to_json(jsonParser &json) const {
    /*
    if(!m_coordinate_space.isIdentity(1e-5))
      json["coordinate_space"]=coordinate_space();
    */
    if(!m_excluded_occs.empty())
      json["excluded_occupants"] = m_excluded_occs;

    DoF::traits(type_name()).to_json(*this, json);
    return json;
  }

}

