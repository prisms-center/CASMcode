#include "casm/clex/ConfigCorrelations.hh"

#include "casm/clex/Clexulator.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/NeighborhoodInfo.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clexulator/ClexParamPack.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

namespace {

/// Return const reference to vector of sequential indices of size >= n
std::vector<unsigned int> const &all_correlation_indices(Index n) {
  static std::vector<unsigned int> all_correlation_indices;
  if (all_correlation_indices.size() < n) {
    all_correlation_indices.reserve(n);
    unsigned int i = all_correlation_indices.size();
    for (; i < n; ++i) {
      all_correlation_indices.push_back(i);
    }
  }
  return all_correlation_indices;
}

}  // namespace

// /// \brief Returns correlations using 'clexulator'. Supercell needs a
// correctly
// /// populated neighbor list.
// Eigen::VectorXd correlations(const ConfigDoF &configdof, const Supercell
// &scel,
//                              Clexulator const &clexulator) {
//   // Size of the supercell will be used for normalizing correlations to a per
//   // primitive cell value
//   int n_unitcells = scel.volume();
//
//   Eigen::VectorXd correlations =
//   Eigen::VectorXd::Zero(clexulator.corr_size());
//
//   // Inform Clexulator of the bitstring
//
//   // Holds contribution to global correlations from a particular Neighborhood
//   Eigen::VectorXd tcorr = correlations;
//
//   for (int unitcell_index = 0; unitcell_index < n_unitcells;
//   unitcell_index++) {
//     // Fill up contributions
//     clexulator.calc_global_corr_contribution(
//         configdof, scel.nlist().sites(unitcell_index).data(),
//         end_ptr(scel.nlist().sites(unitcell_index)), tcorr.data(),
//         end_ptr(tcorr));
//
//     correlations += tcorr;
//   }
//
//   correlations /= (double)n_unitcells;
//
//   return correlations;
// }

/// \brief Returns correlations using 'clexulator'.
Eigen::VectorXd correlations(Configuration const &config,
                             Clexulator const &clexulator) {
  Eigen::VectorXd corr = Eigen::VectorXd::Zero(clexulator.corr_size());
  correlations(corr, config.configdof(), config.supercell().nlist(),
               clexulator);
  return corr;
}

/// \brief Sets correlations using 'clexulator'. Mean of the contribution from
/// every unit cell.
void correlations(Eigen::VectorXd &corr, ConfigDoF const &configdof,
                  SuperNeighborList const &supercell_neighbor_list,
                  Clexulator const &clexulator) {
  auto n = clexulator.corr_size();
  auto const &correlation_indices = all_correlation_indices(n);
  restricted_correlations(corr, configdof, supercell_neighbor_list, clexulator,
                          correlation_indices.data(),
                          correlation_indices.data() + n);
}

/// \brief Sets correlations using 'clexulator', restricted to specified
/// correlation indices. Mean of the contribution from every unit cell.
void restricted_correlations(Eigen::VectorXd &corr, ConfigDoF const &configdof,
                             SuperNeighborList const &supercell_neighbor_list,
                             Clexulator const &clexulator,
                             unsigned int const *corr_indices_begin,
                             unsigned int const *corr_indices_end) {
  restricted_extensive_correlations(corr, configdof, supercell_neighbor_list,
                                    clexulator, corr_indices_begin,
                                    corr_indices_end);
  corr /= (double)supercell_neighbor_list.n_unitcells();
}

/// \brief Returns correlations using 'clexulator'. Sum of the contribution from
/// every unit cell.
void extensive_correlations(Eigen::VectorXd &corr, ConfigDoF const &configdof,
                            SuperNeighborList const &supercell_neighbor_list,
                            Clexulator const &clexulator) {
  auto n = clexulator.corr_size();
  auto const &correlation_indices = all_correlation_indices(n);
  restricted_extensive_correlations(corr, configdof, supercell_neighbor_list,
                                    clexulator, correlation_indices.data(),
                                    correlation_indices.data() + n);
}

/// \brief Sets correlations using 'clexulator', restricted to specified
/// correlation indices. Sum of the contribution from every unit cell.
///
/// \returns Eigen::VectorXd correlations, of size `clexulator.corr_size()`,
/// with zero value for any correlations not in `correlations_indices`.
void restricted_extensive_correlations(
    Eigen::VectorXd &corr, ConfigDoF const &configdof,
    SuperNeighborList const &supercell_neighbor_list,
    Clexulator const &clexulator, unsigned int const *corr_indices_begin,
    unsigned int const *corr_indices_end) {
  int n_corr = clexulator.corr_size();
  int n_unitcells = supercell_neighbor_list.n_unitcells();

  // Holds contribution to global correlations from a particular Neighborhood
  static Eigen::VectorXd tcorr;
  corr.resize(n_corr);
  tcorr.resize(n_corr);

  for (auto it = corr_indices_begin; it != corr_indices_end; ++it) {
    *(corr.data() + *it) = 0.0;
  }

  for (int unitcell_index = 0; unitcell_index < n_unitcells; unitcell_index++) {
    // Fill up contributions
    clexulator.calc_restricted_global_corr_contribution(
        configdof, supercell_neighbor_list.sites(unitcell_index).data(),
        end_ptr(supercell_neighbor_list.sites(unitcell_index)), tcorr.data(),
        end_ptr(tcorr), corr_indices_begin, corr_indices_end);

    for (auto it = corr_indices_begin; it != corr_indices_end; ++it) {
      *(corr.data() + *it) += *(tcorr.data() + *it);
    }
  }
}

/// Returns correlation contribution from a single unit cell, not normalized.
Eigen::VectorXd corr_contribution(Index linear_unitcell_index,
                                  Configuration const &config,
                                  Clexulator const &clexulator) {
  Eigen::VectorXd corr = Eigen::VectorXd::Zero(clexulator.corr_size());
  corr_contribution(corr, linear_unitcell_index, config.configdof(),
                    config.supercell().nlist(), clexulator);
  return corr;
}

/// Returns correlation contribution from a single unit cell, not normalized.
void corr_contribution(Eigen::VectorXd &corr, Index linear_unitcell_index,
                       ConfigDoF const &configdof,
                       SuperNeighborList const &supercell_neighbor_list,
                       Clexulator const &clexulator) {
  int n_corr = clexulator.corr_size();
  int n_unitcells = supercell_neighbor_list.n_unitcells();

  if (linear_unitcell_index >= n_unitcells) {
    std::stringstream msg;
    msg << "Error in point_corr: linear_unitcell_index out of range ("
        << linear_unitcell_index << " >= " << n_unitcells << ")";
    throw std::runtime_error(msg.str());
  }
  corr.resize(n_corr);

  auto const &unitcell_nlist =
      supercell_neighbor_list.sites(linear_unitcell_index);
  clexulator.calc_global_corr_contribution(configdof, unitcell_nlist.data(),
                                           end_ptr(unitcell_nlist), corr.data(),
                                           end_ptr(corr));
}

/// Returns correlation contributions from all unit cells, not normalized.
///
/// \returns Matrix of size (config.supercell().volume(),
/// clexulator.corr_size()). Each row is the contribution of cluster functions
/// associated with a single unit cell.
///
/// Note:
/// - Convert between linear unitcell index (row index) and xtal::UnitCell using
///   SupercellSymInfo::unitcell_index_converter().
/// - Sum over rows is equal to the value returned by extensive_correlations
/// - Average over rows is equal to value returned by correlations
Eigen::MatrixXd all_corr_contribution(Configuration const &config,
                                      Clexulator const &clexulator) {
  int n_unitcells = config.supercell().volume();
  Eigen::MatrixXd corr(n_unitcells, clexulator.corr_size());
  for (Index i = 0; i < n_unitcells; i++) {
    corr.row(i) = corr_contribution(i, config, clexulator);
  }
  return corr;
}

/// \brief Returns point correlations from a single site, normalized by cluster
/// orbit size
Eigen::VectorXd point_corr(Index linear_unitcell_index, Index neighbor_index,
                           const ConfigDoF &configdof, const Supercell &scel,
                           Clexulator const &clexulator) {
  Eigen::VectorXd corr;
  auto n = clexulator.corr_size();
  auto const &correlation_indices = all_correlation_indices(n);
  restricted_point_corr(corr, linear_unitcell_index, neighbor_index, configdof,
                        scel.nlist(), clexulator, correlation_indices.data(),
                        correlation_indices.data() + n);
  return corr;
}

/// \brief Returns point correlations from a single site, normalized by cluster
/// orbit size
Eigen::VectorXd point_corr(Index linear_unitcell_index, Index neighbor_index,
                           const Configuration &config,
                           Clexulator const &clexulator) {
  return point_corr(linear_unitcell_index, neighbor_index, config.configdof(),
                    config.supercell(), clexulator);
}

/// \brief Sets point correlations from a single site, normalized by cluster
/// orbit size, restricted to specified correlations
void restricted_point_corr(Eigen::VectorXd &corr, Index linear_unitcell_index,
                           Index neighbor_index, ConfigDoF const &configdof,
                           SuperNeighborList const &supercell_neighbor_list,
                           Clexulator const &clexulator,
                           unsigned int const *corr_indices_begin,
                           unsigned int const *corr_indices_end) {
  if (linear_unitcell_index >= supercell_neighbor_list.n_unitcells()) {
    std::stringstream msg;
    msg << "Error in restricted_point_corr: linear_unitcell_index out of range "
           "("
        << linear_unitcell_index
        << " >= " << supercell_neighbor_list.n_unitcells() << ")";
    throw std::runtime_error(msg.str());
  }

  int n_corr = clexulator.corr_size();
  corr.resize(n_corr);

  auto const &unitcell_nlist =
      supercell_neighbor_list.sites(linear_unitcell_index);
  int n_point_corr = clexulator.n_point_corr();
  if (neighbor_index >= n_point_corr) {
    std::stringstream msg;
    msg << "Error in restricted_point_corr: neighbor_index out of range ("
        << neighbor_index << " >= " << n_point_corr << ")";
    throw std::runtime_error(msg.str());
  }
  clexulator.calc_restricted_point_corr(
      configdof, unitcell_nlist.data(), end_ptr(unitcell_nlist), neighbor_index,
      corr.data(), end_ptr(corr), corr_indices_begin, corr_indices_end);
}

/// \brief Returns point correlations from all sites, normalized by cluster
/// orbit size
///
/// \returns Matrix of size (n_sites, clexulator.corr_size()), where n_sites =
/// clexulator.n_point_corr() * config.supercell().volume().
///
/// Notes:
/// - The point correlations (each row in the result), described elsewhere in
/// CASM as "flower functions", are the values of all cluster functions that
/// include a particular site, normalized by cluster orbit size.
/// - The interpretation and application differs for periodic vs local cluster
/// expansions
/// - The value clexulator.n_point_corr() is the number of sites for which
/// point correlations can be evaluated (per unit cell), which is the sum of
/// cluster orbit size over all point cluster orbits.
Eigen::MatrixXd all_point_corr(Configuration const &config,
                               Clexulator const &clexulator) {
  Index n_rows = config.supercell().volume() * clexulator.n_point_corr();
  Eigen::MatrixXd corr(n_rows, clexulator.corr_size());
  Index l = 0;
  for (Index j = 0; j < clexulator.n_point_corr(); j++) {
    for (Index i = 0; i < config.supercell().volume(); i++) {
      corr.row(l) = point_corr(i, j, config, clexulator);
      ++l;
    }
  }
  return corr;
}

/// \brief Returns point correlations from all sites, normalized by cluster
/// orbit size, restricted to specified correlations
///
/// \returns Matrix of size (n_sites, clexulator.corr_size()), where n_sites =
/// clexulator.n_point_corr() * config.supercell().volume().
///
/// Notes:
/// - The point correlations (each row in the result), described elsewhere in
/// CASM as "flower functions", are the values of all cluster functions that
/// include a particular site, normalized by cluster orbit size.
/// - The interpretation and application differs for periodic vs local cluster
/// expansions
/// - The value clexulator.n_point_corr() is the number of sites for which
/// point correlations can be evaluated (per unit cell), which is the sum of
/// cluster orbit size over all point cluster orbits.
Eigen::MatrixXd all_restricted_point_corr(
    ConfigDoF const &configdof,
    SuperNeighborList const &supercell_neighbor_list,
    Clexulator const &clexulator, unsigned int const *corr_indices_begin,
    unsigned int const *corr_indices_end) {
  Index n_unitcells = supercell_neighbor_list.n_unitcells();
  Index n_rows = n_unitcells * clexulator.n_point_corr();
  Index n_corr = clexulator.corr_size();
  Eigen::MatrixXd corr = Eigen::MatrixXd::Zero(n_rows, n_corr);
  Eigen::VectorXd tcorr = Eigen::VectorXd::Zero(n_corr);
  Index l = 0;
  for (Index j = 0; j < clexulator.n_point_corr(); j++) {
    for (Index i = 0; i < n_unitcells; i++) {
      restricted_point_corr(tcorr, i, j, configdof, supercell_neighbor_list,
                            clexulator, corr_indices_begin, corr_indices_end);
      corr.row(l) = tcorr;
      ++l;
    }
  }
  return corr;
}

// --- Occupation ---

/// \brief Sets change in (extensive) correlations due to an occupation
/// change, restricted to specified correlations
///
/// \param dcorr, Eigen::VectorXd of change in correlations. Will be set to
/// size `clexulator.corr_size()` if necessary.  Only elements corresponding to
/// indices in `correlations_indices` will be modified.
///
void restricted_delta_corr(Eigen::VectorXd &dcorr, Index linear_site_index,
                           int new_occ, ConfigDoF const &configdof,
                           SuperNeighborList const &supercell_neighbor_list,
                           Clexulator const &clexulator,
                           unsigned int const *corr_indices_begin,
                           unsigned int const *corr_indices_end) {
  int n_corr = clexulator.corr_size();
  dcorr.resize(n_corr);

  Index unitcell_index =
      supercell_neighbor_list.unitcell_index(linear_site_index);
  int neighbor_index =
      supercell_neighbor_list.neighbor_index(linear_site_index);

  auto const &nlist_sites = supercell_neighbor_list.sites(unitcell_index);
  long int const *nlist_begin = nlist_sites.data();
  long int const *nlist_end = end_ptr(nlist_sites);
  double *corr_begin = dcorr.data();
  double *corr_end = end_ptr(dcorr);

  if (!supercell_neighbor_list.overlaps()) {
    int curr_occ = configdof.occ(linear_site_index);
    clexulator.calc_restricted_delta_point_corr(
        configdof, nlist_begin, nlist_end, neighbor_index, curr_occ, new_occ,
        corr_begin, corr_end, corr_indices_begin, corr_indices_end);
  } else {
    static Eigen::VectorXd before;
    before.resize(n_corr);
    Eigen::VectorXd &after = dcorr;

    /// here we have to mutate and unmutate the dof, so in the end it
    /// will not be changed.
    ConfigDoF &mutable_configdof = const_cast<ConfigDoF &>(configdof);

    int curr_occ = configdof.occ(linear_site_index);

    // Calculate before
    clexulator.calc_restricted_point_corr(
        configdof, nlist_begin, nlist_end, neighbor_index, before.data(),
        end_ptr(before), corr_indices_begin, corr_indices_end);

    // Apply change
    mutable_configdof.occ(linear_site_index) = new_occ;

    // Calculate after
    clexulator.calc_restricted_point_corr(
        configdof, nlist_begin, nlist_end, neighbor_index, after.data(),
        end_ptr(after), corr_indices_begin, corr_indices_end);

    // dcorr = after - before
    for (auto it = corr_indices_begin; it != corr_indices_end; ++it) {
      *(dcorr.data() + *it) -= *(before.data() + *it);
    }

    // Unapply changes
    mutable_configdof.occ(linear_site_index) = curr_occ;
  }
}

// --- Local continuous ---

/// \brief Sets change in (extensive) correlations due to a local continuous
/// DoF change, restricted to specified correlations
///
/// \param dcorr, Eigen::VectorXd of change in correlations. Will be set to
/// size `clexulator.corr_size()` if necessary.  Only elements corresponding to
/// indices in `correlations_indices` will be modified.
///
void restricted_delta_corr(Eigen::VectorXd &dcorr, Index linear_site_index,
                           Eigen::VectorXd const &new_value,
                           ConfigDoF const &configdof,
                           SuperNeighborList const &supercell_neighbor_list,
                           LocalContinuousConfigDoFValues const &dof_values,
                           Clexulator const &clexulator,
                           unsigned int const *corr_indices_begin,
                           unsigned int const *corr_indices_end) {
  int n_corr = clexulator.corr_size();
  dcorr.resize(n_corr);

  Index unitcell_index =
      supercell_neighbor_list.unitcell_index(linear_site_index);
  int neighbor_index =
      supercell_neighbor_list.neighbor_index(linear_site_index);

  auto const &nlist_sites = supercell_neighbor_list.sites(unitcell_index);
  long int const *nlist_begin = nlist_sites.data();
  long int const *nlist_end = end_ptr(nlist_sites);

  static Eigen::VectorXd before;
  before.resize(n_corr);
  Eigen::VectorXd &after = dcorr;

  /// here we have to mutate and unmutate the dof, so in the end it
  /// will not be changed.
  LocalContinuousConfigDoFValues &mutable_dof_values =
      const_cast<LocalContinuousConfigDoFValues &>(dof_values);

  Eigen::VectorXd curr_value = dof_values.site_value(linear_site_index);

  // Calculate before
  clexulator.calc_restricted_point_corr(
      configdof, nlist_begin, nlist_end, neighbor_index, before.data(),
      end_ptr(before), corr_indices_begin, corr_indices_end);

  // Apply change
  mutable_dof_values.site_value(linear_site_index) = new_value;

  // Calculate after
  clexulator.calc_restricted_point_corr(
      configdof, nlist_begin, nlist_end, neighbor_index, after.data(),
      end_ptr(after), corr_indices_begin, corr_indices_end);

  // Unapply changes
  mutable_dof_values.site_value(linear_site_index) = curr_value;

  // dcorr = after - before
  for (auto it = corr_indices_begin; it != corr_indices_end; ++it) {
    *(dcorr.data() + *it) -= *(before.data() + *it);
  }
}

// --- Global continuous ---

/// \brief Sets change in (extensive) correlations due to a global continuous
/// DoF change, restricted to specified correlations
///
/// \returns Eigen::VectorXd correlations, of size `clexulator.corr_size()`,
/// with zero value for any correlations not in `correlations_indices`.
void restricted_delta_corr(
    Eigen::VectorXd &dcorr, Eigen::VectorXd const &new_value,
    Eigen::VectorXd const &current_extensive_correlations,
    ConfigDoF const &configdof,
    SuperNeighborList const &supercell_neighbor_list,
    GlobalContinuousConfigDoFValues const &dof_values,
    Clexulator const &clexulator, unsigned int const *corr_indices_begin,
    unsigned int const *corr_indices_end) {
  /// here we have to mutate and unmutate the dof, so in the end it
  /// will not be changed.
  GlobalContinuousConfigDoFValues &mutable_dof_values =
      const_cast<GlobalContinuousConfigDoFValues &>(dof_values);

  Eigen::VectorXd curr_value = dof_values.values();

  // Apply change
  mutable_dof_values.set_values(new_value);

  // Calculate after (set directly to dcorr)
  restricted_extensive_correlations(dcorr, configdof, supercell_neighbor_list,
                                    clexulator, corr_indices_begin,
                                    corr_indices_end);

  // Unapply changes
  mutable_dof_values.set_values(curr_value);

  Eigen::VectorXd const &before = current_extensive_correlations;
  // dcorr = after - before
  for (auto it = corr_indices_begin; it != corr_indices_end; ++it) {
    *(dcorr.data() + *it) -= *(before.data() + *it);
  }
}

/// Return xtal::UnitCellCoord for each row in `all_point_corr`
std::vector<xtal::UnitCellCoord> make_all_point_corr_unitcellcoord(
    NeighborhoodInfo const &neighborhood_info,
    SupercellSymInfo const &sym_info) {
  auto T = sym_info.transformation_matrix_to_super();
  auto f = sym_info.unitcellcoord_index_converter();
  auto lattice_points = xtal::make_lattice_points(T);
  auto const &origin_unitcellcoords =
      neighborhood_info.point_corr_unitcellcoord;

  std::vector<xtal::UnitCellCoord> result;
  for (xtal::UnitCellCoord unitcellcoord : origin_unitcellcoords) {
    for (xtal::UnitCell trans : lattice_points) {
      result.push_back(f.bring_within(unitcellcoord + trans));
    }
  }
  return result;
}

/// Return coordinate (xtal::Coordinate) for each row in `all_point_corr`
///
/// - Coordinates are referenced to the supercell lattice
std::vector<xtal::Coordinate> make_all_point_corr_coordinates(
    xtal::BasicStructure const &prim, NeighborhoodInfo const &neighborhood_info,
    SupercellSymInfo const &sym_info) {
  xtal::Lattice const &supercell_lattice = sym_info.supercell_lattice();
  std::vector<xtal::UnitCellCoord> uccoords =
      make_all_point_corr_unitcellcoord(neighborhood_info, sym_info);

  std::vector<xtal::Coordinate> result;
  for (xtal::UnitCellCoord const &uccoord : uccoords) {
    xtal::Coordinate coord = uccoord.coordinate(prim);
    coord.set_lattice(supercell_lattice, CART);
    result.push_back(coord);
  }
  return result;
}

/// Return Cartesian coordinate (as rows) for each row in `all_point_corr`
Eigen::MatrixXd make_all_point_corr_cart_coordinates(
    xtal::BasicStructure const &prim, NeighborhoodInfo const &neighborhood_info,
    SupercellSymInfo const &sym_info) {
  std::vector<xtal::Coordinate> coords =
      make_all_point_corr_coordinates(prim, neighborhood_info, sym_info);
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(coords.size(), 3);
  Index l = 0;
  for (xtal::Coordinate const &coord : coords) {
    result.row(l) = coord.const_cart();
    l++;
  }
  return result;
}

/// Return fractional coordinate (as rows) for each row in `all_point_corr`
///
/// - Fractional coordinates are in terms of the supercell lattice vectors
Eigen::MatrixXd make_all_point_corr_frac_coordinates(
    xtal::BasicStructure const &prim, NeighborhoodInfo const &neighborhood_info,
    SupercellSymInfo const &sym_info) {
  std::vector<xtal::Coordinate> coords =
      make_all_point_corr_coordinates(prim, neighborhood_info, sym_info);
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(coords.size(), 3);
  Index l = 0;
  for (xtal::Coordinate coord : coords) {
    result.row(l) = coord.const_frac();
    l++;
  }
  return result;
}

/// Return asymmetric unit index for each row in `all_point_corr`
///
/// Note:
/// - Asymmetric unit index is in range [0, n_point_orbits), where
///   n_point_orbits is the number of point orbits in the cluster orbits used to
///   construct the clexulator
std::vector<Index> make_all_point_corr_asymmetric_unit_indices(
    NeighborhoodInfo const &neighborhood_info,
    SupercellSymInfo const &sym_info) {
  std::vector<Index> result;
  Index volume = sym_info.unitcell_index_converter().total_sites();
  auto const &indices = neighborhood_info.point_corr_asymmetric_unit_indices;
  for (Index j = 0; j < indices.size(); ++j) {
    for (Index i = 0; i < volume; ++i) {
      result.push_back(indices[j]);
    }
  }
  return result;
}

/// Return linear unitcell index for each row in `all_point_corr`
std::vector<Index> make_all_point_corr_linear_unitcell_indices(
    NeighborhoodInfo const &neighborhood_info,
    SupercellSymInfo const &sym_info) {
  std::vector<Index> result;
  Index volume = sym_info.unitcell_index_converter().total_sites();
  for (Index j = 0; j < neighborhood_info.n_point_corr; ++j) {
    for (Index i = 0; i < volume; ++i) {
      result.push_back(i);
    }
  }
  return result;
}

/// \brief Returns gradient correlations using 'clexulator', with respect to DoF
/// 'dof_type'
Eigen::MatrixXd gradcorrelations(ConfigDoF const &configdof,
                                 Supercell const &scel,
                                 Clexulator const &clexulator, DoFKey &key) {
  clexulator::ClexParamKey paramkey;
  clexulator::ClexParamKey corr_key(clexulator.param_pack().key("corr"));
  clexulator::ClexParamKey dof_key;
  if (key == "occ") {
    paramkey = clexulator.param_pack().key("diff/corr/" + key + "_site_func");
    dof_key = clexulator.param_pack().key("occ_site_func");
  } else {
    paramkey = clexulator.param_pack().key("diff/corr/" + key + "_var");
    dof_key = clexulator.param_pack().key(key + "_var");
  }

  std::string em_corr, em_dof;
  em_corr = clexulator.param_pack().eval_mode(corr_key);
  em_dof = clexulator.param_pack().eval_mode(dof_key);

  // this const_cast is not great...
  // but it seems like the only place passing const Clexulator is a problem and
  // it is not actually changing clexulator before/after this function
  const_cast<Clexulator &>(clexulator)
      .param_pack()
      .set_eval_mode(corr_key, "DIFF");
  const_cast<Clexulator &>(clexulator)
      .param_pack()
      .set_eval_mode(dof_key, "DIFF");

  Eigen::MatrixXd gcorr;
  Index scel_vol = scel.volume();
  if (DoF::BasicTraits(key).global()) {
    Eigen::MatrixXd gcorr_func = configdof.global_dof(key).values();
    gcorr.setZero(gcorr_func.size(), clexulator.corr_size());
    // Holds contribution to global correlations from a particular Neighborhood

    // std::vector<double> corr(clexulator.corr_size(), 0.0);
    for (int v = 0; v < scel_vol; v++) {
      // Fill up contributions
      clexulator.calc_global_corr_contribution(configdof,
                                               scel.nlist().sites(v).data(),
                                               end_ptr(scel.nlist().sites(v)));

      for (Index c = 0; c < clexulator.corr_size(); ++c)
        gcorr.col(c) += clexulator.param_pack().read(paramkey(c));
    }
  } else {
    Eigen::MatrixXd gcorr_func;
    gcorr.setZero(configdof.local_dof(key).values().size(),
                  clexulator.corr_size());
    // Holds contribution to global correlations from a particular Neighborhood
    Index l;
    for (int v = 0; v < scel_vol; v++) {
      // Fill up contributions
      clexulator.calc_global_corr_contribution(configdof,
                                               scel.nlist().sites(v).data(),
                                               end_ptr(scel.nlist().sites(v)));

      for (Index c = 0; c < clexulator.corr_size(); ++c) {
        gcorr_func = clexulator.param_pack().read(paramkey(c));

        for (Index n = 0; n < clexulator.nlist_size(); ++n) {
          l = scel.nlist().sites(v)[n];
          // for(Index i=0; i<gcorr_func.cols(); ++i){
          gcorr.block(l * gcorr_func.rows(), c, gcorr_func.rows(), 1) +=
              gcorr_func.col(n);
          // std::cout << "Block: (" << l * gcorr_func.rows() << ", " << c << ",
          // " << gcorr_func.rows() << ", " << 1 << ") += " <<
          // gcorr_func.col(n).transpose() << "\n";
          //}
        }
      }
    }
  }
  const_cast<Clexulator &>(clexulator)
      .param_pack()
      .set_eval_mode(corr_key, em_corr);
  const_cast<Clexulator &>(clexulator)
      .param_pack()
      .set_eval_mode(dof_key, em_dof);

  return gcorr;
}

/// \brief Returns gradient correlations using 'clexulator', with respect to DoF
/// 'dof_type'
Eigen::MatrixXd gradcorrelations(const Configuration &config,
                                 Clexulator const &clexulator, DoFKey &key) {
  return gradcorrelations(config.configdof(), config.supercell(), clexulator,
                          key);
}

}  // namespace CASM
