#include "casm/basis_set/DoF.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {

  void DoF::register_traits(BasicTraits const &_type) {
    _traits_map().insert(_type);
    //std::cout << "Adding DoFType " << _type.name() << "\nMap now holds:\n";
    //for(auto  const & d : _traits_map()){
    //std::cout << "  " << d.name() << "\n";
    //}
    //std::cout << "\n";
  }

  //********************************************************************

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




}

