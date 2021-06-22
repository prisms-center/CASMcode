#include "casm/clex/ConfigCorrelations.hh"

#include "casm/clex/ClexParamPack.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

/// \brief Returns correlations using 'clexulator'. Supercell needs a correctly
/// populated neighbor list.
Eigen::VectorXd correlations(const ConfigDoF &configdof, const Supercell &scel,
                             Clexulator const &clexulator) {
  // Size of the supercell will be used for normalizing correlations to a per
  // primitive cell value
  int n_unitcells = scel.volume();

  Eigen::VectorXd correlations = Eigen::VectorXd::Zero(clexulator.corr_size());

  // Inform Clexulator of the bitstring

  // Holds contribution to global correlations from a particular neighborhood
  Eigen::VectorXd tcorr = correlations;
  // std::vector<double> corr(clexulator.corr_size(), 0.0);

  for (int unitcell_index = 0; unitcell_index < n_unitcells; unitcell_index++) {
    // Fill up contributions
    clexulator.calc_global_corr_contribution(
        configdof, scel.nlist().sites(unitcell_index).data(),
        end_ptr(scel.nlist().sites(unitcell_index)), tcorr.data(),
        end_ptr(tcorr));

    correlations += tcorr;
  }

  correlations /= (double)n_unitcells;

  return correlations;
}

/// \brief Returns correlations using 'clexulator', restricted to specified
/// correlation indices. Supercell needs a correctly populated neighbor list.
///
/// \returns Eigen::VectorXd correlations, of size `clexulator.corr_size()`,
/// with zero value for any correlations not in `correlations_indices`.
Eigen::VectorXd restricted_correlations(
    const ConfigDoF &configdof, const Supercell &scel,
    Clexulator const &clexulator,
    std::vector<unsigned int> const &correlation_indices) {
  // Size of the supercell will be used for normalizing correlations to a per
  // primitive cell value
  int n_unitcells = scel.volume();

  Eigen::VectorXd correlations = Eigen::VectorXd::Zero(clexulator.corr_size());

  // Holds contribution to global correlations from a particular neighborhood
  Eigen::VectorXd tcorr = correlations;

  auto corr_indices_begin = correlation_indices.data();
  auto corr_indices_end = corr_indices_begin + correlation_indices.size();

  for (int unitcell_index = 0; unitcell_index < n_unitcells; unitcell_index++) {
    // Fill up contributions
    clexulator.calc_restricted_global_corr_contribution(
        configdof, scel.nlist().sites(unitcell_index).data(),
        end_ptr(scel.nlist().sites(unitcell_index)), tcorr.data(),
        end_ptr(tcorr), corr_indices_begin, corr_indices_end);

    correlations += tcorr;
  }

  correlations /= (double)n_unitcells;

  return correlations;
}

/// \brief Returns correlations using 'clexulator'.
Eigen::VectorXd correlations(const Configuration &config,
                             Clexulator const &clexulator) {
  return correlations(config.configdof(), config.supercell(), clexulator);
}

/// \brief Returns correlations using 'clexulator', restricted to specified
/// correlation indices.
///
/// \returns Eigen::VectorXd correlations, of size `clexulator.corr_size()`,
/// with zero value for any correlations not in `correlations_indices`.
Eigen::VectorXd restricted_correlations(
    const Configuration &config, Clexulator const &clexulator,
    std::vector<unsigned int> const &correlation_indices) {
  return restricted_correlations(config.configdof(), config.supercell(),
                                 clexulator, correlation_indices);
}

/// \brief Returns correlations using 'clexulator'. Sum of the contribution from
/// every unit cell.
Eigen::VectorXd extensive_correlations(Configuration const &config,
                                       Clexulator const &clexulator) {
  ConfigDoF const &configdof = config.configdof();
  SuperNeighborList const &supercell_neighbor_list = config.supercell().nlist();

  // number of unit cells
  int n_unitcells = config.supercell().volume();

  Eigen::VectorXd correlations = Eigen::VectorXd::Zero(clexulator.corr_size());

  // Holds contribution to global correlations from a particular neighborhood
  Eigen::VectorXd tcorr = correlations;

  for (int unitcell_index = 0; unitcell_index < n_unitcells; unitcell_index++) {
    // Fill up contributions
    clexulator.calc_global_corr_contribution(
        configdof, supercell_neighbor_list.sites(unitcell_index).data(),
        end_ptr(supercell_neighbor_list.sites(unitcell_index)), tcorr.data(),
        end_ptr(tcorr));
    correlations += tcorr;
  }
  return correlations;
}

/// \brief Returns correlations using 'clexulator', restricted to specified
/// correlation indices. Sum of the contribution from every unit cell.
///
/// \returns Eigen::VectorXd correlations, of size `clexulator.corr_size()`,
/// with zero value for any correlations not in `correlations_indices`.
Eigen::VectorXd restricted_extensive_correlations(
    Configuration const &config, Clexulator const &clexulator,
    std::vector<unsigned int> const &correlation_indices) {
  ConfigDoF const &configdof = config.configdof();
  SuperNeighborList const &supercell_neighbor_list = config.supercell().nlist();

  // number of unit cells
  int n_unitcells = config.supercell().volume();

  Eigen::VectorXd correlations = Eigen::VectorXd::Zero(clexulator.corr_size());

  // Holds contribution to global correlations from a particular neighborhood
  Eigen::VectorXd tcorr = correlations;

  auto corr_indices_begin = correlation_indices.data();
  auto corr_indices_end = corr_indices_begin + correlation_indices.size();

  for (int unitcell_index = 0; unitcell_index < n_unitcells; unitcell_index++) {
    // Fill up contributions
    clexulator.calc_restricted_global_corr_contribution(
        configdof, supercell_neighbor_list.sites(unitcell_index).data(),
        end_ptr(supercell_neighbor_list.sites(unitcell_index)), tcorr.data(),
        end_ptr(tcorr), corr_indices_begin, corr_indices_end);
    correlations += tcorr;
  }
  return correlations;
}

/// Returns correlation contribution from a single unit cell, not normalized.
///
/// Supercell needs a correctly populated neighbor list.
Eigen::VectorXd corr_contribution(Index linear_unitcell_index,
                                  const ConfigDoF &configdof,
                                  const Supercell &scel,
                                  Clexulator const &clexulator) {
  if (linear_unitcell_index >= scel.volume()) {
    std::stringstream msg;
    msg << "Error in point_corr: linear_unitcell_index out of range ("
        << linear_unitcell_index << " >= " << scel.volume() << ")";
    throw std::runtime_error(msg.str());
  }
  Eigen::VectorXd correlations = Eigen::VectorXd::Zero(clexulator.corr_size());

  auto const &unitcell_nlist = scel.nlist().sites(linear_unitcell_index);
  clexulator.calc_global_corr_contribution(
      configdof, unitcell_nlist.data(), end_ptr(unitcell_nlist),
      correlations.data(), end_ptr(correlations));

  return correlations;
}

/// \brief Returns correlations using 'clexulator'.
Eigen::VectorXd corr_contribution(Index linear_unitcell_index,
                                  const Configuration &config,
                                  Clexulator const &clexulator) {
  return corr_contribution(linear_unitcell_index, config.configdof(),
                           config.supercell(), clexulator);
}

/// Returns correlation contributions from all unit cells, not normalized.
Eigen::MatrixXd all_corr_contribution(const Configuration &config,
                                      Clexulator const &clexulator) {
  Eigen::MatrixXd corr(config.supercell().volume(), clexulator.corr_size());
  for (Index i = 0; i < config.supercell().volume(); i++) {
    corr.row(i) = corr_contribution(i, config, clexulator);
  }
  return corr;
}

/// \brief Returns point correlations from a single site, normalized by cluster
/// orbit size
Eigen::VectorXd point_corr(Index linear_unitcell_index, Index neighbor_index,
                           const ConfigDoF &configdof, const Supercell &scel,
                           Clexulator const &clexulator) {
  if (linear_unitcell_index >= scel.volume()) {
    std::stringstream msg;
    msg << "Error in point_corr: linear_unitcell_index out of range ("
        << linear_unitcell_index << " >= " << scel.volume() << ")";
    throw std::runtime_error(msg.str());
  }

  Eigen::VectorXd correlations = Eigen::VectorXd::Zero(clexulator.corr_size());

  auto const &unitcell_nlist = scel.nlist().sites(linear_unitcell_index);
  if (neighbor_index >= clexulator.n_point_corr()) {
    std::stringstream msg;
    msg << "Error in point_corr: neighbor_index out of range ("
        << neighbor_index << " >= " << clexulator.n_point_corr() << ")";
    throw std::runtime_error(msg.str());
  }
  clexulator.calc_point_corr(configdof, unitcell_nlist.data(),
                             end_ptr(unitcell_nlist), neighbor_index,
                             correlations.data(), end_ptr(correlations));

  return correlations;
}

/// \brief Returns point correlations from a single site, normalized by cluster
/// orbit size
Eigen::VectorXd point_corr(Index linear_unitcell_index, Index neighbor_index,
                           const Configuration &config,
                           Clexulator const &clexulator) {
  return point_corr(linear_unitcell_index, neighbor_index, config.configdof(),
                    config.supercell(), clexulator);
}

/// \brief Returns point correlations from all sites, normalized by cluster
/// orbit size
Eigen::MatrixXd all_point_corr(const Configuration &config,
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

// --- Occupation ---

/// \brief Returns change in (extensive) correlations due to an occupation
/// change
Eigen::VectorXd delta_corr(Index linear_site_index, int new_occ,
                           Configuration const &configuration,
                           Clexulator const &clexulator) {
  Eigen::VectorXd dcorr{Eigen::VectorXd::Zero(clexulator.corr_size())};

  ConfigDoF const &configdof = configuration.configdof();
  SuperNeighborList const &supercell_neighbor_list =
      configuration.supercell().nlist();

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
    clexulator.calc_delta_point_corr(configdof, nlist_begin, nlist_end,
                                     neighbor_index, curr_occ, new_occ,
                                     corr_begin, corr_end);
  } else {
    Eigen::VectorXd before{Eigen::VectorXd::Zero(dcorr.size())};
    Eigen::VectorXd after{Eigen::VectorXd::Zero(dcorr.size())};

    /// here we have to mutate and unmutate the dof, so in the end it
    /// will not be changed.
    ConfigDoF &mutable_configdof = const_cast<ConfigDoF &>(configdof);

    int curr_occ = configdof.occ(linear_site_index);

    // Calculate before
    clexulator.calc_point_corr(configdof, nlist_begin, nlist_end,
                               neighbor_index, before.data(), end_ptr(before));

    // Apply change
    mutable_configdof.occ(linear_site_index) = new_occ;

    // Calculate after
    clexulator.calc_point_corr(configdof, nlist_begin, nlist_end,
                               neighbor_index, after.data(), end_ptr(after));

    dcorr = after - before;

    // Unapply changes
    mutable_configdof.occ(linear_site_index) = curr_occ;
  }

  return dcorr;
}

/// \brief Returns change in (extensive) correlations due to an occupation
/// change, restricted to specified correlations
///
/// \returns Eigen::VectorXd correlations, of size `clexulator.corr_size()`,
/// with zero value for any correlations not in `correlations_indices`.
Eigen::VectorXd restricted_delta_corr(
    Index linear_site_index, int new_occ, Configuration const &configuration,
    Clexulator const &clexulator,
    std::vector<unsigned int> const &correlation_indices) {
  Eigen::VectorXd dcorr{Eigen::VectorXd::Zero(clexulator.corr_size())};
  restricted_delta_corr(dcorr, linear_site_index, new_occ, configuration,
                        clexulator, correlation_indices);
  return dcorr;
}

/// \brief Returns change in (extensive) correlations due to an occupation
/// change, restricted to specified correlations
///
/// \returns Eigen::VectorXd correlations, of size `clexulator.corr_size()`,
/// with values set only for any correlations in `correlations_indices`.
void restricted_delta_corr(
    Eigen::VectorXd &dcorr, Index linear_site_index, int new_occ,
    Configuration const &configuration, Clexulator const &clexulator,
    std::vector<unsigned int> const &correlation_indices) {
  if (dcorr.size() != clexulator.corr_size()) {
    dcorr = Eigen::VectorXd::Zero(clexulator.corr_size());
  }

  ConfigDoF const &configdof = configuration.configdof();
  SuperNeighborList const &supercell_neighbor_list =
      configuration.supercell().nlist();

  Index unitcell_index =
      supercell_neighbor_list.unitcell_index(linear_site_index);
  int neighbor_index =
      supercell_neighbor_list.neighbor_index(linear_site_index);

  auto const &nlist_sites = supercell_neighbor_list.sites(unitcell_index);
  long int const *nlist_begin = nlist_sites.data();
  long int const *nlist_end = end_ptr(nlist_sites);
  double *corr_begin = dcorr.data();
  double *corr_end = end_ptr(dcorr);
  unsigned int const *corr_index_begin = correlation_indices.data();
  unsigned int const *corr_index_end = end_ptr(correlation_indices);

  if (!supercell_neighbor_list.overlaps()) {
    int curr_occ = configuration.occ(linear_site_index);
    clexulator.calc_restricted_delta_point_corr(
        configdof, nlist_begin, nlist_end, neighbor_index, curr_occ, new_occ,
        corr_begin, corr_end, corr_index_begin, corr_index_end);
  } else {
    Eigen::VectorXd before{Eigen::VectorXd::Zero(dcorr.size())};
    Eigen::VectorXd after{Eigen::VectorXd::Zero(dcorr.size())};

    /// here we have to mutate and unmutate the dof, so in the end it
    /// will not be changed.
    ConfigDoF &mutable_configdof = const_cast<ConfigDoF &>(configdof);

    int curr_occ = configdof.occ(linear_site_index);

    // Calculate before
    clexulator.calc_restricted_point_corr(
        configdof, nlist_begin, nlist_end, neighbor_index, before.data(),
        end_ptr(before), corr_index_begin, corr_index_end);

    // Apply change
    mutable_configdof.occ(linear_site_index) = new_occ;

    // Calculate after
    clexulator.calc_restricted_point_corr(
        configdof, nlist_begin, nlist_end, neighbor_index, after.data(),
        end_ptr(after), corr_index_begin, corr_index_end);

    dcorr = after - before;

    // Unapply changes
    mutable_configdof.occ(linear_site_index) = curr_occ;
  }
}

// --- Local continuous ---

/// \brief Returns change in (extensive) correlations due to a local continuous
/// DoF change
Eigen::VectorXd delta_corr(Index linear_site_index,
                           Eigen::VectorXd const &new_value, DoFKey const &key,
                           Configuration const &configuration,
                           Clexulator const &clexulator) {
  Eigen::VectorXd dcorr{Eigen::VectorXd::Zero(clexulator.corr_size())};

  ConfigDoF const &configdof = configuration.configdof();
  SuperNeighborList const &supercell_neighbor_list =
      configuration.supercell().nlist();
  LocalContinuousConfigDoFValues const &dof_values = configdof.local_dof(key);

  Index unitcell_index =
      supercell_neighbor_list.unitcell_index(linear_site_index);
  int neighbor_index =
      supercell_neighbor_list.neighbor_index(linear_site_index);

  auto const &nlist_sites = supercell_neighbor_list.sites(unitcell_index);
  long int const *nlist_begin = nlist_sites.data();
  long int const *nlist_end = end_ptr(nlist_sites);

  Eigen::VectorXd before{Eigen::VectorXd::Zero(dcorr.size())};
  Eigen::VectorXd after{Eigen::VectorXd::Zero(dcorr.size())};

  /// here we have to mutate and unmutate the dof, so in the end it
  /// will not be changed.
  LocalContinuousConfigDoFValues &mutable_dof_values =
      const_cast<LocalContinuousConfigDoFValues &>(dof_values);

  Eigen::VectorXd curr_value = dof_values.site_value(linear_site_index);

  // Calculate before
  clexulator.calc_point_corr(configdof, nlist_begin, nlist_end, neighbor_index,
                             before.data(), end_ptr(before));

  // Apply change
  mutable_dof_values.site_value(linear_site_index) = new_value;

  // Calculate after
  clexulator.calc_point_corr(configdof, nlist_begin, nlist_end, neighbor_index,
                             after.data(), end_ptr(after));

  dcorr = after - before;

  // Unapply changes
  mutable_dof_values.site_value(linear_site_index) = curr_value;

  return dcorr;
}

/// \brief Returns change in (extensive) correlations due to a local continuous
/// DoF change, restricted to specified correlations
///
/// \returns Eigen::VectorXd correlations, of size `clexulator.corr_size()`,
/// with zero value for any correlations not in `correlations_indices`.
Eigen::VectorXd restricted_delta_corr(
    Index linear_site_index, Eigen::VectorXd const &new_value,
    DoFKey const &key, Configuration const &configuration,
    Clexulator const &clexulator,
    std::vector<unsigned int> const &correlation_indices) {
  Eigen::VectorXd dcorr{Eigen::VectorXd::Zero(clexulator.corr_size())};

  ConfigDoF const &configdof = configuration.configdof();
  SuperNeighborList const &supercell_neighbor_list =
      configuration.supercell().nlist();
  LocalContinuousConfigDoFValues const &dof_values = configdof.local_dof(key);

  Index unitcell_index =
      supercell_neighbor_list.unitcell_index(linear_site_index);
  int neighbor_index =
      supercell_neighbor_list.neighbor_index(linear_site_index);

  auto const &nlist_sites = supercell_neighbor_list.sites(unitcell_index);
  long int const *nlist_begin = nlist_sites.data();
  long int const *nlist_end = end_ptr(nlist_sites);
  unsigned int const *corr_index_begin = correlation_indices.data();
  unsigned int const *corr_index_end = end_ptr(correlation_indices);

  Eigen::VectorXd before{Eigen::VectorXd::Zero(dcorr.size())};
  Eigen::VectorXd after{Eigen::VectorXd::Zero(dcorr.size())};

  /// here we have to mutate and unmutate the dof, so in the end it
  /// will not be changed.
  LocalContinuousConfigDoFValues &mutable_dof_values =
      const_cast<LocalContinuousConfigDoFValues &>(dof_values);

  Eigen::VectorXd curr_value = dof_values.site_value(linear_site_index);

  // Calculate before
  clexulator.calc_restricted_point_corr(
      configdof, nlist_begin, nlist_end, neighbor_index, before.data(),
      end_ptr(before), corr_index_begin, corr_index_end);

  // Apply change
  mutable_dof_values.site_value(linear_site_index) = new_value;

  // Calculate after
  clexulator.calc_restricted_point_corr(
      configdof, nlist_begin, nlist_end, neighbor_index, after.data(),
      end_ptr(after), corr_index_begin, corr_index_end);

  dcorr = after - before;

  // Unapply changes
  mutable_dof_values.site_value(linear_site_index) = curr_value;

  return dcorr;
}

// --- Global continuous ---

/// \brief Returns change in (extensive) correlations due to a global continuous
/// DoF change
Eigen::VectorXd delta_corr(
    Eigen::VectorXd const &new_value, DoFKey const &key,
    Eigen::VectorXd const &current_extensive_correlations,
    Configuration const &configuration, Clexulator const &clexulator) {
  GlobalContinuousConfigDoFValues const &dof_values =
      configuration.configdof().global_dof(key);

  /// here we have to mutate and unmutate the dof, so in the end it
  /// will not be changed.
  GlobalContinuousConfigDoFValues &mutable_dof_values =
      const_cast<GlobalContinuousConfigDoFValues &>(dof_values);

  Eigen::VectorXd curr_value = dof_values.values();

  // Apply change
  mutable_dof_values.set_values(new_value);

  // Calculate after
  Eigen::VectorXd after = extensive_correlations(configuration, clexulator);

  // Unapply changes
  mutable_dof_values.set_values(curr_value);

  return after - current_extensive_correlations;
}

/// \brief Returns change in (extensive) correlations due to a global continuous
/// DoF change, restricted to specified correlations
///
/// \returns Eigen::VectorXd correlations, of size `clexulator.corr_size()`,
/// with zero value for any correlations not in `correlations_indices`.
Eigen::VectorXd restricted_delta_corr(
    Eigen::VectorXd const &new_value, DoFKey const &key,
    Eigen::VectorXd const &current_extensive_correlations,
    Configuration const &configuration, Clexulator const &clexulator,
    std::vector<unsigned int> const &correlation_indices) {
  GlobalContinuousConfigDoFValues const &dof_values =
      configuration.configdof().global_dof(key);

  /// here we have to mutate and unmutate the dof, so in the end it
  /// will not be changed.
  GlobalContinuousConfigDoFValues &mutable_dof_values =
      const_cast<GlobalContinuousConfigDoFValues &>(dof_values);

  Eigen::VectorXd curr_value = dof_values.values();

  // Apply change
  mutable_dof_values.set_values(new_value);

  // Calculate after
  Eigen::VectorXd after = restricted_extensive_correlations(
      configuration, clexulator, correlation_indices);

  // Unapply changes
  mutable_dof_values.set_values(curr_value);

  return after - current_extensive_correlations;
}

/// Return a vector of xtal::Coordinate for each row in `all_point_corr`
std::vector<xtal::UnitCellCoord> make_all_point_corr_unitcellcoord(
    Index n_point_corr, xtal::BasicStructure const &prim,
    SupercellSymInfo const &sym_info, SuperNeighborList const &scel_nlist) {
  // from the input:
  auto f = sym_info.unitcellcoord_index_converter();
  Index volume = sym_info.superlattice().size();

  std::vector<xtal::UnitCellCoord> result;
  for (Index neighbor_index = 0; neighbor_index < n_point_corr;
       neighbor_index++) {
    for (Index unitcell_index = 0; unitcell_index < volume; unitcell_index++) {
      auto const &neighbor_list = scel_nlist.sites(unitcell_index);
      Index linear_site_index = neighbor_list[neighbor_index];
      result.push_back(f(linear_site_index));
    }
  }
  return result;
}

/// Return coordinates (xtal::Coordinate) for each row in `all_point_corr`
///
/// - Coordinates are referenced to the supercell lattice
std::vector<xtal::Coordinate> make_all_point_corr_coordinates(
    Index n_point_corr, xtal::BasicStructure const &prim,
    SupercellSymInfo const &sym_info, SuperNeighborList const &scel_nlist) {
  std::vector<xtal::UnitCellCoord> uccoords = make_all_point_corr_unitcellcoord(
      n_point_corr, prim, sym_info, scel_nlist);

  xtal::Lattice const &supercell_lattice = sym_info.supercell_lattice();

  std::vector<xtal::Coordinate> result;
  for (xtal::UnitCellCoord const &uccoord : uccoords) {
    xtal::Coordinate coord = uccoord.coordinate(prim);
    coord.set_lattice(supercell_lattice, CART);
    result.push_back(coord);
  }
  return result;
}

/// Return Cartesian coordinates (as rows) for each row in `all_point_corr`
Eigen::MatrixXd make_all_point_corr_cart_coordinates(
    Index n_point_corr, xtal::BasicStructure const &prim,
    SupercellSymInfo const &sym_info, SuperNeighborList const &scel_nlist) {
  std::vector<xtal::Coordinate> coords =
      make_all_point_corr_coordinates(n_point_corr, prim, sym_info, scel_nlist);
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(coords.size(), 3);
  Index l = 0;
  for (xtal::Coordinate const &coord : coords) {
    result.row(l) = coord.const_cart();
    l++;
  }
  return result;
}

/// Return fractional coordinates (as rows) for each row in `all_point_corr`
///
/// - Fractional coordinates are in terms of the supercell lattice vectors
Eigen::MatrixXd make_all_point_corr_frac_coordinates(
    Index n_point_corr, xtal::BasicStructure const &prim,
    SupercellSymInfo const &sym_info, SuperNeighborList const &scel_nlist) {
  std::vector<xtal::Coordinate> coords =
      make_all_point_corr_coordinates(n_point_corr, prim, sym_info, scel_nlist);
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(coords.size(), 3);
  Index l = 0;
  for (xtal::Coordinate coord : coords) {
    result.row(l) = coord.const_frac();
    l++;
  }
  return result;
}

/// \brief Returns gradient correlations using 'clexulator', with respect to DoF
/// 'dof_type'
Eigen::MatrixXd gradcorrelations(const ConfigDoF &configdof,
                                 const Supercell &scel,
                                 Clexulator const &clexulator, DoFKey &key) {
  ClexParamKey paramkey;
  ClexParamKey corr_key(clexulator.param_pack().key("corr"));
  ClexParamKey dof_key;
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
    // Holds contribution to global correlations from a particular neighborhood

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
    // Holds contribution to global correlations from a particular neighborhood
    Index l;
    for (int v = 0; v < scel_vol; v++) {
      // Fill up contributions
      clexulator.calc_global_corr_contribution(configdof,
                                               scel.nlist().sites(v).data(),
                                               end_ptr(scel.nlist().sites(v)));

      for (Index c = 0; c < clexulator.corr_size(); ++c) {
        gcorr_func = clexulator.param_pack().read(paramkey(c));

        for (Index n = 0; n < scel.nlist().sites(v).size(); ++n) {
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
