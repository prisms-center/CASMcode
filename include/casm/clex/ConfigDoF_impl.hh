#ifndef CASM_clex_ConfigDoF_impl
#define CASM_clex_ConfigDoF_impl

#include <vector>

#include "casm/clex/ConfigDoF.hh"

namespace CASM {

/// ConfigDoF constructor
///
/// \param _N_sublat Number of sublattices in corresponding prim
/// \param _N_vol Supercell volume, as multiples of the prim volume
/// \param global_dof_info GlobalInfoContainerType is an iterable container of
///        value_type std::pair<DoFKey,ContinuousDoFInfo>.
/// \param local_dof_info LocalInfoContainerType is an iterable container of
///        value_type std::pair<DoFKey,std::vector<ContinuousDoFInfo>  >
//
/// Note: Typically Structure has already been included when a ConfigDoF is
/// constructed, in which case it is best to use `make_configdof` from
/// clex/ConfigDoFTools.hh.
///
template <typename GlobalInfoContainerType, typename LocalInfoContainerType>
ConfigDoF::ConfigDoF(Index _N_sublat, Index _N_vol,
                     GlobalInfoContainerType const &global_dof_info,
                     LocalInfoContainerType const &local_dof_info,
                     std::vector<SymGroupRepID> const &occ_symrep_IDs,
                     double _tol)
    : m_occupation(DoF::BasicTraits("occ"), _N_sublat, _N_vol, occ_symrep_IDs),
      m_tol(_tol) {
  for (auto const &dof : global_dof_info) {
    DoF::BasicTraits ttraits(dof.first);

    if (!ttraits.global())
      throw std::runtime_error(
          "Attempting to initialize ConfigDoF global value using local DoF " +
          dof.first);
    m_global_dofs.emplace(
        dof.first, GlobalContinuousConfigDoFValues(ttraits, _N_sublat, _N_vol,
                                                   dof.second));
  }
  for (auto const &dof : local_dof_info) {
    DoF::BasicTraits ttraits(dof.first);
    if (_N_sublat == 0) continue;
    if (ttraits.global())
      throw std::runtime_error(
          "Attempting to initialize ConfigDoF local value using global DoF " +
          dof.first);
    if (_N_sublat != dof.second.size()) {
      throw std::runtime_error(
          "Attempting to initialize ConfigDoF local value '" + dof.first +
          "' with improperly initialized parameter 'local_dof_info'.");
    }

    m_local_dofs.emplace(
        dof.first,
        LocalContinuousConfigDoFValues(ttraits, _N_sublat, _N_vol, dof.second));
  }
}

}  // namespace CASM

#endif
