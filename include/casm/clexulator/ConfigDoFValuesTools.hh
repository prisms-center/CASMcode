#ifndef CASM_clexulator_ConfigDoFValuesTools
#define CASM_clexulator_ConfigDoFValuesTools
#include <algorithm>

#include "casm/clexulator/ConfigDoFValues.hh"

namespace CASM {
namespace clexulator {

// Note: DoFSetType template parameter requires: dim(), basis(), inv_basis()

/// Returns the block of DoF values from one sublattice
template <typename Derived>
Eigen::Block<Derived> sublattice_block(Eigen::MatrixBase<Derived> &M,
                                       Index sublattice_index, Index N_volume);

/// Returns the block of DoF values from one sublattice
template <typename Derived>
const Eigen::Block<const Derived> sublattice_block(
    Eigen::MatrixBase<Derived> const &M, Index sublattice_index,
    Index N_volume);

/// \brief Determine local continuous DoF matrix row dimension
///
/// For local continuous DoF values stored in a matrix, the expected matrix
/// shape is (matrix_dim x number of sites), where matrix_dim is the maximum
/// of the individual sublattice DoFSet dimensions.
///
/// \param dof_info Vector of DoFSetType, which must have a `DoFSetType::dim()`
/// member
template <typename DoFSetType>
Index max_dim(std::vector<DoFSetType> const &dof_info);

/// Convert local DoF values from prim basis to standard basis
template <typename DoFSetType>
Eigen::MatrixXd local_to_standard_values(
    Eigen::MatrixXd const &dof_values, Index N_sublat, Index N_volume,
    std::vector<DoFSetType> const &dof_info);

/// Convert local DoF values from standard basis to prim basis
template <typename DoFSetType>
Eigen::MatrixXd local_from_standard_values(
    Eigen::MatrixXd const &standard_values, Index N_sublat, Index N_volume,
    std::vector<DoFSetType> const &dof_info);

/// Convert global DoF values from prim basis to standard basis
template <typename DoFSetType>
Eigen::VectorXd global_to_standard_values(Eigen::VectorXd const &dof_values,
                                          DoFSetType const &dof_info);

/// Convert global DoF values from standard basis to prim basis
template <typename DoFSetType>
Eigen::VectorXd global_from_standard_values(
    Eigen::VectorXd const &standard_values, DoFSetType const &dof_info);

/// \brief Make ConfigDoFValues, expressed in the standard basis, initialized
/// with value zero
template <typename GlobalDoFSetType, typename LocalDoFSetType>
ConfigDoFValues make_default_standard_config_dof_values(
    Index N_sublat, Index N_volume,
    std::map<DoFKey, GlobalDoFSetType> const &global_dof_info,
    std::map<DoFKey, std::vector<LocalDoFSetType>> const &local_dof_info);

/// \brief Make ConfigDoFValues, expressed in the prim basis, initialized
/// with value zero
template <typename GlobalDoFSetType, typename LocalDoFSetType>
ConfigDoFValues make_default_config_dof_values(
    Index N_sublat, Index N_volume,
    std::map<DoFKey, GlobalDoFSetType> const &global_dof_info,
    std::map<DoFKey, std::vector<LocalDoFSetType>> const &local_dof_info);

/// Convert ConfigDoFValues expressed in prim basis to standard basis
template <typename GlobalDoFSetType, typename LocalDoFSetType>
ConfigDoFValues to_standard_values(
    ConfigDoFValues const &dof_values, Index N_sublat, Index N_volume,
    std::map<DoFKey, GlobalDoFSetType> const &global_dof_info,
    std::map<DoFKey, std::vector<LocalDoFSetType>> const &local_dof_info);

/// Convert ConfigDoFValues expressed in standard basis to prim basis
template <typename GlobalDoFSetType, typename LocalDoFSetType>
ConfigDoFValues from_standard_values(
    ConfigDoFValues const &standard_values, Index N_sublat, Index N_volume,
    std::map<DoFKey, GlobalDoFSetType> const &global_dof_info,
    std::map<DoFKey, std::vector<LocalDoFSetType>> const &local_dof_info);

// --- Inline & template definitions ---

/// Returns the block of DoF values from one sublattice
template <typename Derived>
Eigen::Block<Derived> sublattice_block(Eigen::MatrixBase<Derived> &M,
                                       Index sublattice_index, Index N_volume) {
  return Eigen::Block<Derived>(M.derived(), 0, sublattice_index * N_volume,
                               M.rows(), N_volume);
}

/// Returns the block of DoF values from one sublattice
template <typename Derived>
const Eigen::Block<Derived const> sublattice_block(
    Eigen::MatrixBase<Derived> const &M, Index sublattice_index,
    Index N_volume) {
  return Eigen::Block<Derived const>(
      M.derived(), 0, sublattice_index * N_volume, M.rows(), N_volume);
}

/// \brief Determine local continuous DoF matrix row dimension
///
/// For local continuous DoF values stored in a matrix, the expected matrix
/// shape is (matrix_dim x number of sites), where matrix_dim is the maximum
/// of the individual sublattice DoFSet dimensions.
///
/// \param dof_info Vector of DoFSetType, which must have a `DoFSetType::dim()`
/// member
template <typename DoFSetType>
Index max_dim(std::vector<DoFSetType> const &dof_info) {
  Index dim = 0;
  for (auto const &info : dof_info) {
    dim = std::max(dim, info.dim());
  }
  return dim;
}

/// Convert local DoF values from prim basis to standard basis
template <typename DoFSetType>
Eigen::MatrixXd local_to_standard_values(
    Eigen::MatrixXd const &dof_values, Index N_sublat, Index N_volume,
    std::vector<DoFSetType> const &dof_info) {
  Index rows = dof_info.front().basis().rows();
  Eigen::MatrixXd standard_values(rows, dof_values.cols());
  for (Index b = 0; b < N_sublat; ++b) {
    standard_values.block(0, b * N_volume, rows, N_volume) =
        dof_info[b].basis() *
        sublattice_block(dof_values, b, N_volume).topRows(dof_info[b].dim());
  }
  return standard_values;
}

/// Convert local DoF values from standard basis to prim basis
template <typename DoFSetType>
Eigen::MatrixXd local_from_standard_values(
    Eigen::MatrixXd const &standard_values, Index N_sublat, Index N_volume,
    std::vector<DoFSetType> const &dof_info) {
  Index N_sites = N_volume * N_sublat;
  if (standard_values.rows() != dof_info.front().basis().rows() ||
      standard_values.cols() != N_sites) {
    std::stringstream msg;
    msg << "Invalid standard values input size in "
           "local_dof_values_from_standard_basis: "
        << "Expected rows=" << dof_info.front().basis().rows()
        << ", received rows=" << standard_values.rows()
        << ", expected cols=" << N_sites
        << ", received cols=" << standard_values.cols();
    throw std::runtime_error(msg.str());
  }
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(max_dim(dof_info), N_sites);
  for (Index b = 0; b < N_sublat; ++b)
    sublattice_block(result, b, N_volume).topRows(dof_info[b].dim()) =
        dof_info[b].inv_basis() *
        sublattice_block(standard_values, b, N_volume);
  return result;
}

/// Convert global DoF values from prim basis to standard basis
template <typename DoFSetType>
Eigen::VectorXd global_to_standard_values(Eigen::VectorXd const &dof_values,
                                          DoFSetType const &dof_info) {
  return dof_info.basis() * dof_values;
}

/// Convert global DoF values from standard basis to prim basis
template <typename DoFSetType>
Eigen::VectorXd global_from_standard_values(
    Eigen::VectorXd const &standard_values, DoFSetType const &dof_info) {
  if (standard_values.size() != dof_info.basis().rows()) {
    std::stringstream msg;
    msg << "Invalid standard values input size in "
           "global_dof_values_from_standard_basis: "
        << "Expected size=" << dof_info.basis().rows()
        << ", received size=" << standard_values.size();
    throw std::runtime_error(msg.str());
  }
  return dof_info.inv_basis() * standard_values;
}

/// \brief Make ConfigDoFValues, expressed in the standard basis, initialized
/// with value zero
template <typename GlobalDoFSetType, typename LocalDoFSetType>
ConfigDoFValues make_default_standard_config_dof_values(
    Index N_sublat, Index N_volume,
    std::map<DoFKey, GlobalDoFSetType> const &global_dof_info,
    std::map<DoFKey, std::vector<LocalDoFSetType>> const &local_dof_info) {
  Index N_sites = N_volume * N_sublat;
  ConfigDoFValues standard_values;
  standard_values.occupation = Eigen::VectorXi::Zero(N_sites);

  for (auto const &pair : local_dof_info) {
    std::string const &dof_type = pair.first;
    auto const &dof_info = pair.second;
    Index standard_dim = dof_info.front().basis().rows();
    standard_values.local_dof_values.emplace(
        dof_type, Eigen::MatrixXd::Zero(standard_dim, N_sites));
  }

  for (auto const &pair : global_dof_info) {
    std::string const &dof_type = pair.first;
    auto const &dof_info = pair.second;
    Index standard_dim = dof_info.basis().rows();
    standard_values.local_dof_values.emplace(
        dof_type, Eigen::VectorXd::Zero(standard_dim));
  }

  return standard_values;
}

/// \brief Make ConfigDoFValues, expressed in the prim basis, initialized
/// with value zero
template <typename GlobalDoFSetType, typename LocalDoFSetType>
ConfigDoFValues make_default_config_dof_values(
    Index N_sublat, Index N_volume,
    std::map<DoFKey, GlobalDoFSetType> const &global_dof_info,
    std::map<DoFKey, std::vector<LocalDoFSetType>> const &local_dof_info) {
  Index N_sites = N_volume * N_sublat;
  ConfigDoFValues dof_values;
  dof_values.occupation = Eigen::VectorXi::Zero(N_sites);

  for (auto const &pair : local_dof_info) {
    std::string const &dof_type = pair.first;
    auto const &dof_info = pair.second;
    dof_values.local_dof_values.emplace(
        dof_type, Eigen::MatrixXd::Zero(max_dim(dof_info), N_sites));
  }

  for (auto const &pair : global_dof_info) {
    std::string const &dof_type = pair.first;
    auto const &dof_info = pair.second;
    dof_values.local_dof_values.emplace(dof_type,
                                        Eigen::VectorXd::Zero(dof_info.dim()));
  }

  return dof_values;
}

/// Convert ConfigDoFValues expressed in prim basis to standard basis
template <typename GlobalDoFSetType, typename LocalDoFSetType>
ConfigDoFValues to_standard_values(
    ConfigDoFValues const &dof_values, Index N_sublat, Index N_volume,
    std::map<DoFKey, GlobalDoFSetType> const &global_dof_info,
    std::map<DoFKey, std::vector<LocalDoFSetType>> const &local_dof_info) {
  ConfigDoFValues standard_values;
  standard_values.occupation = dof_values.occupation;

  for (auto const &pair : dof_values.local_dof_values) {
    std::string const &dof_type = pair.first;
    Eigen::MatrixXd const &values = pair.second;
    auto it = local_dof_info.find(dof_type);
    if (it == local_dof_info.end()) {
      std::stringstream msg;
      msg << "Error in clexulator::to_standard_values: Incompatible local DoF";
      throw std::runtime_error(msg.str());
    }
    auto const &dof_info = it->second;
    standard_values.local_dof_values.emplace(
        dof_type,
        local_to_standard_values(values, N_sublat, N_volume, dof_info));
  }

  for (auto const &pair : dof_values.global_dof_values) {
    std::string const &dof_type = pair.first;
    Eigen::VectorXd const &values = pair.second;
    auto it = global_dof_info.find(dof_type);
    if (it == global_dof_info.end()) {
      std::stringstream msg;
      msg << "Error in clexulator::to_standard_values: Incompatible local DoF";
      throw std::runtime_error(msg.str());
    }
    auto const &dof_info = it->second;
    standard_values.local_dof_values.emplace(
        dof_type, global_to_standard_values(values, dof_info));
  }

  return standard_values;
}

/// Convert ConfigDoFValues expressed in standard basis to prim basis
template <typename GlobalDoFSetType, typename LocalDoFSetType>
ConfigDoFValues from_standard_values(
    ConfigDoFValues const &standard_values, Index N_sublat, Index N_volume,
    std::map<DoFKey, GlobalDoFSetType> const &global_dof_info,
    std::map<DoFKey, std::vector<LocalDoFSetType>> const &local_dof_info) {
  ConfigDoFValues dof_values;
  dof_values.occupation = dof_values.occupation;

  for (auto const &pair : standard_values.local_dof_values) {
    std::string const &dof_type = pair.first;
    Eigen::MatrixXd const &values = pair.second;
    auto it = local_dof_info.find(dof_type);
    if (it == local_dof_info.end()) {
      std::stringstream msg;
      msg << "Error in clexulator::to_standard_values: Incompatible local DoF";
      throw std::runtime_error(msg.str());
    }
    auto const &dof_info = it->second;
    dof_values.local_dof_values.emplace(
        dof_type,
        local_from_standard_values(values, N_sublat, N_volume, dof_info));
  }

  for (auto const &pair : standard_values.global_dof_values) {
    std::string const &dof_type = pair.first;
    Eigen::VectorXd const &values = pair.second;
    auto it = global_dof_info.find(dof_type);
    if (it == global_dof_info.end()) {
      std::stringstream msg;
      msg << "Error in clexulator::to_standard_values: Incompatible local DoF";
      throw std::runtime_error(msg.str());
    }
    auto const &dof_info = it->second;
    dof_values.local_dof_values.emplace(
        dof_type, global_from_standard_values(values, dof_info));
  }

  return dof_values;
}

}  // namespace clexulator
}  // namespace CASM

#endif
