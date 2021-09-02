#ifndef CASM_clex_ConfigDoF_impl
#define CASM_clex_ConfigDoF_impl

#include <vector>

#include "casm/clex/ConfigDoF.hh"
#include "casm/clexulator/ConfigDoFValuesTools.hh"

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
    : m_N_sublat(_N_sublat),
      m_N_vol(_N_vol),
      m_dof_values(),
      m_occ_symrep_IDs(occ_symrep_IDs),
      m_tol(_tol) {
  Index N_sites = m_N_sublat * m_N_vol;
  m_dof_values.occupation = Eigen::VectorXi::Zero(N_sites);

  for (auto const &dof : global_dof_info) {
    std::string const &dof_name = dof.first;
    DoFSetInfo const &dof_info = dof.second;
    DoF::BasicTraits ttraits(dof_name);

    if (!ttraits.global())
      throw std::runtime_error(
          "Attempting to initialize ConfigDoF global value using local DoF " +
          dof.first);
    auto result = m_dof_values.global_dof_values.emplace(
        dof_name, Eigen::VectorXd::Zero(dof_info.dim()));
    Eigen::VectorXd &values = result.first->second;
    m_global_dofs.try_emplace(dof_name, ttraits, _N_sublat, _N_vol, dof_info,
                              values);
  }
  for (auto const &dof : local_dof_info) {
    std::string const &dof_name = dof.first;
    std::vector<DoFSetInfo> const &dof_info = dof.second;
    DoF::BasicTraits ttraits(dof_name);
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

    auto result = m_dof_values.local_dof_values.emplace(
        dof_name,
        Eigen::MatrixXd::Zero(clexulator::max_dim(dof_info), N_sites));
    Eigen::MatrixXd &values = result.first->second;
    m_local_dofs.try_emplace(dof_name, ttraits, _N_sublat, _N_vol, dof_info,
                             values);
  }
}

}  // namespace CASM

#endif
