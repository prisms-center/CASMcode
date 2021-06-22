#ifndef CASM_ConfigCorrelations
#define CASM_ConfigCorrelations

#include "casm/crystallography/DoFDecl.hh"
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {

namespace xtal {
class BasicStructure;
class Coordinate;
class UnitCellCoord;
}  // namespace xtal

class Clexulator;
class ConfigDoF;
class Configuration;
class SuperNeighborList;
class Supercell;
class SupercellSymInfo;

// --- Correlations ---

/// \brief Returns correlations using 'clexulator'. Supercell needs a correctly
/// populated neighbor list.
Eigen::VectorXd correlations(const ConfigDoF &configdof, const Supercell &scel,
                             Clexulator const &clexulator);

/// \brief Returns correlations using 'clexulator', restricted to specified
/// correlation indices.
Eigen::VectorXd restricted_correlations(
    const ConfigDoF &configdof, const Supercell &scel,
    Clexulator const &clexulator,
    std::vector<unsigned int> const &correlation_indices);

/// \brief Returns correlations using 'clexulator'.
Eigen::VectorXd correlations(const Configuration &config,
                             Clexulator const &clexulator);

/// \brief Returns correlations using 'clexulator', restricted to specified
/// correlation indices.
Eigen::VectorXd restricted_correlations(
    const Configuration &config, Clexulator const &clexulator,
    std::vector<unsigned int> const &correlation_indices);

/// \brief Returns correlations using 'clexulator'. Sum of the contribution from
/// every unit cell.
Eigen::VectorXd extensive_correlations(Configuration const &config,
                                       Clexulator const &clexulator);

/// \brief Returns correlations using 'clexulator', restricted to specified
/// correlation indices. Sum of the contribution from every unit cell.
Eigen::VectorXd restricted_extensive_correlations(
    Configuration const &config, Clexulator const &clexulator,
    std::vector<unsigned int> const &correlation_indices);

// --- Correlations, contribution of a particular unit cell ---

/// Returns correlation contribution from a single unit cell, not normalized.
Eigen::VectorXd corr_contribution(Index linear_unitcell_index,
                                  const ConfigDoF &configdof,
                                  const Supercell &scel,
                                  Clexulator const &clexulator);

/// \brief Returns correlations using 'clexulator'.
Eigen::VectorXd corr_contribution(Index linear_unitcell_index,
                                  const Configuration &config,
                                  Clexulator const &clexulator);

/// Returns correlation contributions from all unit cells, not normalized.
Eigen::MatrixXd all_corr_contribution(const Configuration &config,
                                      Clexulator const &clexulator);

// --- Point correlations ---

/// \brief Returns point correlations from a single site, normalized by cluster
/// orbit size
Eigen::VectorXd point_corr(Index linear_unitcell_index, Index neighbor_index,
                           const ConfigDoF &configdof, const Supercell &scel,
                           Clexulator const &clexulator);

/// \brief Returns point correlations from a single site, normalized by cluster
/// orbit size
Eigen::VectorXd point_corr(Index linear_unitcell_index, Index neighbor_index,
                           const Configuration &config,
                           Clexulator const &clexulator);

/// \brief Returns point correlations from all sites, normalized by cluster
/// orbit size
Eigen::MatrixXd all_point_corr(const Configuration &config,
                               Clexulator const &clexulator);

// --- Occupation ---

/// \brief Returns change in (extensive) correlations due to an occupation
/// change
Eigen::VectorXd delta_corr(Index linear_site_index, int new_occ,
                           Configuration const &configuration,
                           Clexulator const &clexulator);

/// \brief Returns change in (extensive) correlations due to an occupation
/// change, restricted to specified correlations
Eigen::VectorXd restricted_delta_corr(
    Index linear_site_index, int new_occ, Configuration const &configuration,
    Clexulator const &clexulator,
    std::vector<unsigned int> const &correlation_indices);

// --- Local continuous ---

/// \brief Returns change in (extensive) correlations due to a local continuous
/// DoF change
Eigen::VectorXd delta_corr(Index linear_site_index,
                           Eigen::VectorXd const &new_value, DoFKey const &key,
                           Configuration const &configuration,
                           Clexulator const &clexulator);

/// \brief Returns change in (extensive) correlations due to a local continuous
/// DoF change, restricted to specified correlations
Eigen::VectorXd restricted_delta_corr(
    Index linear_site_index, Eigen::VectorXd const &new_value,
    DoFKey const &key, Configuration const &configuration,
    Clexulator const &clexulator,
    std::vector<unsigned int> const &correlation_indices);

// --- Global continuous ---

/// \brief Returns change in (extensive) correlations due to a global
/// continuous DoF change
Eigen::VectorXd delta_corr(
    Eigen::VectorXd const &new_value, DoFKey const &key,
    Eigen::VectorXd const &current_extensive_correlations,
    Configuration const &configuration, Clexulator const &clexulator);

/// \brief Returns change in (extensive) correlations due to a global
/// continuous DoF change, restricted to specified correlations
Eigen::VectorXd restricted_delta_corr(
    Eigen::VectorXd const &new_value, DoFKey const &key,
    Eigen::VectorXd const &current_extensive_correlations,
    Configuration const &configuration, Clexulator const &clexulator,
    std::vector<unsigned int> const &correlation_indices);

// --- Coordinates for sites in `all_point_corr` ---

/// Return a vector of xtal::Coordinate for each row in `all_point_corr`
std::vector<xtal::UnitCellCoord> make_all_point_corr_unitcellcoord(
    Index n_point_corr, xtal::BasicStructure const &prim,
    SupercellSymInfo const &sym_info, SuperNeighborList const &scel_nlist);

/// Return coordinates (xtal::Coordinate) for each row in `all_point_corr`
std::vector<xtal::Coordinate> make_all_point_corr_coordinates(
    Index n_point_corr, xtal::BasicStructure const &prim,
    SupercellSymInfo const &sym_info, SuperNeighborList const &scel_nlist);

/// Return Cartesian coordinates (as rows) for each row in `all_point_corr`
Eigen::MatrixXd make_all_point_corr_cart_coordinates(
    Index n_point_corr, xtal::BasicStructure const &prim,
    SupercellSymInfo const &sym_info, SuperNeighborList const &scel_nlist);

/// Return fractional coordinates (as rows) for each row in `all_point_corr`
Eigen::MatrixXd make_all_point_corr_frac_coordinates(
    Index n_point_corr, xtal::BasicStructure const &prim,
    SupercellSymInfo const &sym_info, SuperNeighborList const &scel_nlist);

// --- Gradient correlations ---

/// \brief Returns gradient correlations using 'clexulator', with respect to DoF
/// 'dof_type'
Eigen::MatrixXd gradcorrelations(const ConfigDoF &configdof,
                                 const Supercell &scel,
                                 Clexulator const &clexulator, DoFKey &key);

/// \brief Returns gradient correlations using 'clexulator', with respect to DoF
/// 'dof_type'
Eigen::MatrixXd gradcorrelations(const Configuration &config,
                                 Clexulator const &clexulator, DoFKey &key);

}  // namespace CASM

#endif
