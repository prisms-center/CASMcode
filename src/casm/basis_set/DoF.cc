#include "casm/basis_set/DoF.hh"
#include "casm/symmetry/SymGroup.hh"

namespace CASM {
  namespace DoFType {
    namespace Local {
      static std::map<std::string, BasicTraits> &_traits_map() {
        static std::map<std::string, BasicTraits> _static_traits_map;
        return _static_traits_map;
      }
    }

    //********************************************************************

    void register_traits(BasicTraits const &_type) {
      Local::_traits_map().emplace(_type.name(), _type);
      //std::cout << "Adding DoFType " << _type.name() << "\nMap now holds:\n";
      //for(auto  const & d : _traits_map()){
      //std::cout << "  " << d.name() << "\n";
      //}
      //std::cout << "\n";
    }

    //********************************************************************

    BasicTraits const &traits(std::string const &_type_name) {
      auto it = Local::_traits_map().find(_type_name);
      if(it == Local::_traits_map().end()) {
        throw std::runtime_error("Could not find DoF Traits for DoF type" + _type_name);
      }
      return *it;
    }

  }

  namespace DoF_impl {
    //********************************************************************

    Base::Base(Base::BasicTraits const &_traits,
               std::string const &_var_name,
               Index _ID) :
      m_type_name(_traits.type_name()),
      m_var_name(_var_name),
      m_dof_ID(_ID),
      m_ID_lock(false) {

      DoFType::register_traits(_traits);

    }
  }

  //********************************************************************
  void DiscreteDoF::allocate_symrep(SymGroup const &_group) const {
    if(!m_symrep_ID.empty() && !m_symrep_ID.is_identity())
      throw std::runtime_error("In DiscreteDoF::allocate_symrep(), representation has already been allocated for this DoF.");

    set_symrep_ID(_group.allocate_representation());
  }



}

