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

}

