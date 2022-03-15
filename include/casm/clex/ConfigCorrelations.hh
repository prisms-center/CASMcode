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

namespace clexulator {
class SuperNeighborList;
}
using clexulator::SuperNeighborList;

class Clexulator;
class ConfigDoF;
class Configuration;
class GlobalContinuousConfigDoFValues;
class LocalContinuousConfigDoFValues;
struct NeighborhoodInfo;
class Supercell;
class SupercellSymInfo;

// --- Correlations ---

/// \brief Returns correlations using 'clexulator'.
Eigen::VectorXd correlations(Configuration const &config,
                             Clexulator const &clexulator);

/// \brief Sets correlations using 'clexulator'. Mean of the contribution from
/// every unit cell.
void correlations(Eigen::VectorXd &corr, ConfigDoF const &configdof,
                  SuperNeighborList const &supercell_neighbor_list,
                  Clexulator const &clexulator);

/// \brief Sets correlations using 'clexulator', restricted to specified
/// correlation indices. Mean of the contribution from every unit cell.
void restricted_correlations(Eigen::VectorXd &corr, ConfigDoF const &configdof,
                             SuperNeighborList const &supercell_neighbor_list,
                             Clexulator const &clexulator,
                             unsigned int const *corr_indices_begin,
                             unsigned int const *corr_indices_end);

/// \brief Sets correlations using 'clexulator'. Sum of the contribution from
/// every unit cell.
void extensive_correlations(Eigen::VectorXd &corr, ConfigDoF const &configdof,
                            SuperNeighborList const &supercell_neighbor_list,
                            Clexulator const &clexulator);

/// \brief Sets correlations using 'clexulator', restricted to specified
/// correlation indices. Sum of the contribution from every unit cell.
void restricted_extensive_correlations(
    Eigen::VectorXd &corr, ConfigDoF const &configdof,
    SuperNeighborList const &supercell_neighbor_list,
    Clexulator const &clexulator, unsigned int const *corr_indices_begin,
    unsigned int const *corr_indices_end);

// --- Correlations, contribution of a particular unit cell ---

/// Returns correlation contribution from a single unit cell, not normalized.
Eigen::VectorXd corr_contribution(Index linear_unitcell_index,
                                  Configuration const &config,
                                  Clexulator const &clexulator);

/// Returns correlation contribution from a single unit cell, not normalized.
void corr_contribution(Eigen::VectorXd &corr, Index linear_unitcell_index,
                       ConfigDoF const &configdof,
                       SuperNeighborList const &supercell_neighbor_list,
                       Clexulator const &clexulator);

/// Returns correlation contributions from all unit cells, not normalized.
Eigen::MatrixXd all_corr_contribution(Configuration const &config,
                                      Clexulator const &clexulator);

// --- Point correlations ---

/// \brief Returns point correlations from a single site, normalized by cluster
/// orbit size
Eigen::VectorXd point_corr(Index linear_unitcell_index, Index neighbor_index,
                           ConfigDoF const &configdof, Supercell const &scel,
                           Clexulator const &clexulator);

/// \brief Returns point correlations from a single site, normalized by cluster
/// orbit size
Eigen::VectorXd point_corr(Index linear_unitcell_index, Index neighbor_index,
                           Configuration const &config,
                           Clexulator const &clexulator);

/// \brief Sets point correlations from a single site, normalized by cluster
/// orbit size, restricted to specified correlations
void restricted_point_corr(Eigen::VectorXd &corr, Index linear_unitcell_index,
                           Index neighbor_index, ConfigDoF const &configdof,
                           SuperNeighborList const &supercell_neighbor_list,
                           Clexulator const &clexulator,
                           unsigned int const *corr_indices_begin,
                           unsigned int const *corr_indices_end);

/// \brief Returns point correlations from all sites, normalized by cluster
/// orbit size
Eigen::MatrixXd all_point_corr(Configuration const &config,
                               Clexulator const &clexulator);

/// \brief Returns point correlations from all sites, normalized by cluster
/// orbit size, restricted to specified correlations
Eigen::MatrixXd all_restricted_point_corr(
    ConfigDoF const &configdof,
    SuperNeighborList const &supercell_neighbor_list,
    Clexulator const &clexulator, unsigned int const *corr_indices_begin,
    unsigned int const *corr_indices_end);

// --- Occupation ---

/// \brief Sets change in (extensive) correlations due to an occupation
/// change, restricted to specified correlations
void restricted_delta_corr(Eigen::VectorXd &dcorr, Index linear_site_index,
                           int new_occ, ConfigDoF const &configdof,
                           SuperNeighborList const &supercell_neighbor_list,
                           Clexulator const &clexulator,
                           unsigned int const *corr_indices_begin,
                           unsigned int const *corr_indices_end);

// --- Local continuous ---

/// \brief Sets change in (extensive) correlations due to a local continuous
/// DoF change, restricted to specified correlations
void restricted_delta_corr(Eigen::VectorXd &dcorr, Index linear_site_index,
                           Eigen::VectorXd const &new_value,
                           ConfigDoF const &configdof,
                           SuperNeighborList const &supercell_neighbor_list,
                           LocalContinuousConfigDoFValues const &dof_values,
                           Clexulator const &clexulator,
                           unsigned int const *corr_indices_begin,
                           unsigned int const *corr_indices_end);

// --- Global continuous ---

/// \brief Sets change in (extensive) correlations due to a global
/// continuous DoF change, restricted to specified correlations
void restricted_delta_corr(
    Eigen::VectorXd &dcorr, Eigen::VectorXd const &new_value,
    Eigen::VectorXd const &current_extensive_correlations,
    ConfigDoF const &configdof,
    SuperNeighborList const &supercell_neighbor_list,
    GlobalContinuousConfigDoFValues const &dof_values,
    Clexulator const &clexulator, unsigned int const *corr_indices_begin,
    unsigned int const *corr_indices_end);

// --- Coordinates for sites in `all_point_corr` ---

/// Return xtal::UnitCellCoord for each row in `all_point_corr`
std::vector<xtal::UnitCellCoord> make_all_point_corr_unitcellcoord(
    NeighborhoodInfo const &Neighborhood_info,
    SupercellSymInfo const &sym_info);

/// Return coordinate (xtal::Coordinate) for each row in `all_point_corr`
std::vector<xtal::Coordinate> make_all_point_corr_coordinates(
    xtal::BasicStructure const &prim, NeighborhoodInfo const &Neighborhood_info,
    SupercellSymInfo const &sym_info);

/// Return Cartesian coordinate (as rows) for each row in `all_point_corr`
Eigen::MatrixXd make_all_point_corr_cart_coordinates(
    xtal::BasicStructure const &prim, NeighborhoodInfo const &Neighborhood_info,
    SupercellSymInfo const &sym_info);

/// Return fractional coordinate (as rows) for each row in `all_point_corr`
Eigen::MatrixXd make_all_point_corr_frac_coordinates(
    xtal::BasicStructure const &prim, NeighborhoodInfo const &Neighborhood_info,
    SupercellSymInfo const &sym_info);

/// Return asymmetric unit index for each row in `all_point_corr`
std::vector<Index> make_all_point_corr_asymmetric_unit_indices(
    NeighborhoodInfo const &Neighborhood_info,
    SupercellSymInfo const &sym_info);

/// Return linear unitcell index for each row in `all_point_corr`
std::vector<Index> make_all_point_corr_linear_unitcell_indices(
    NeighborhoodInfo const &Neighborhood_info,
    SupercellSymInfo const &sym_info);

// --- Gradient correlations ---

/// \brief Returns gradient correlations using 'clexulator', with respect to DoF
/// 'dof_type'
Eigen::MatrixXd gradcorrelations(ConfigDoF const &configdof,
                                 Supercell const &scel,
                                 Clexulator const &clexulator, DoFKey &key);

/// \brief Returns gradient correlations using 'clexulator', with respect to DoF
/// 'dof_type'
Eigen::MatrixXd gradcorrelations(Configuration const &config,
                                 Clexulator const &clexulator, DoFKey &key);

}  // namespace CASM

#endif
