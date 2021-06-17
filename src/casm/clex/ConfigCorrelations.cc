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
  int scel_vol = scel.volume();

  Eigen::VectorXd correlations = Eigen::VectorXd::Zero(clexulator.corr_size());

  // Inform Clexulator of the bitstring

  // Holds contribution to global correlations from a particular neighborhood
  Eigen::VectorXd tcorr = correlations;
  // std::vector<double> corr(clexulator.corr_size(), 0.0);

  for (int v = 0; v < scel_vol; v++) {
    // Fill up contributions
    clexulator.calc_global_corr_contribution(
        configdof, scel.nlist().sites(v).data(), end_ptr(scel.nlist().sites(v)),
        tcorr.data(), end_ptr(tcorr));

    correlations += tcorr;
  }

  correlations /= (double)scel_vol;

  return correlations;
}

/// \brief Returns correlations using 'clexulator'.
Eigen::VectorXd correlations(const Configuration &config,
                             Clexulator const &clexulator) {
  return correlations(config.configdof(), config.supercell(), clexulator);
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
