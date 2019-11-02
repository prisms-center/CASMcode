#include "casm/basis_set/DoF.hh"
#include "casm/symmetry/SymGroup.hh"

namespace CASM {
  namespace DoF {
    //********************************************************************

    Base::Base(BasicTraits const &_traits,
               std::string const &_var_name,
               Index _ID) :
      m_traits(_traits),
      m_var_name(_var_name),
      m_dof_ID(_ID),
      m_ID_lock(false) {

      //DoFType::register_traits(_traits);

    }
  }

  //********************************************************************
  void DiscreteDoF::allocate_symrep(SymGroup const &_group) const {
    if(!m_symrep_ID.empty() && !m_symrep_ID.is_identity())
      throw std::runtime_error("In DiscreteDoF::allocate_symrep(), representation has already been allocated for this DoF.");

    set_symrep_ID(_group.allocate_representation());
  }



}

