#include "casm/clex/ConfigIO.hh"

#include <functional>

#include "casm/app/ClexDescription.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/dataformatter/DataFormatterTools_impl.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/Calculable.hh"
#include "casm/clex/ClexBasisSpecs.hh"
#include "casm/clex/ConfigCorrelations.hh"
#include "casm/clex/ConfigIOHull.hh"
#include "casm/clex/ConfigIOLocalCorr.hh"
#include "casm/clex/ConfigIONovelty.hh"
#include "casm/clex/ConfigIOStrain.hh"
#include "casm/clex/ConfigIOStrucScore.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/clex/MappedProperties.hh"
#include "casm/clex/NeighborhoodInfo.hh"
#include "casm/clex/Norm.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/clex/io/json/Configuration_json_io.hh"
#include "casm/clex/io/stream/Configuration_stream_io.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/io/SimpleStructureIO.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"
#include "casm/crystallography/io/VaspIO.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/database/PropertiesDatabase.hh"
#include "casm/database/Selected.hh"

namespace CASM {

template class BaseDatumFormatter<Configuration>;
template class DataFormatterOperator<bool, std::string, Configuration>;
template class DataFormatterOperator<bool, bool, Configuration>;
template class DataFormatterOperator<bool, double, Configuration>;
template class DataFormatterOperator<double, double, Configuration>;
template class DataFormatterOperator<Index, double, Configuration>;
template class DataFormatter<Configuration>;
template bool DataFormatter<Configuration>::evaluate_as_scalar<bool>(
    Configuration const &) const;
template double DataFormatter<Configuration>::evaluate_as_scalar<double>(
    Configuration const &) const;
template class DataFormatterDictionary<Configuration>;

namespace ConfigIO_impl {

/// \brief Expects arguments of the form 'name' or 'name(Au)', 'name(Pt)', etc.
bool MolDependent::parse_args(const std::string &args) {
  if (args.size() > 0) m_mol_names.push_back(args);
  return true;
}

/// \brief Adds index rules corresponding to the parsed args
bool MolDependent::init(const Configuration &_tmplt) const {
  auto struc_mol = xtal::struc_molecule_name(_tmplt.primclex().prim());

  if (m_mol_names.size() == 0) {
    for (Index i = 0; i < struc_mol.size(); i++) {
      _add_rule(std::vector<Index>({i}));
      m_mol_names.push_back(struc_mol[i]);
    }
  } else {
    for (Index n = 0; n < m_mol_names.size(); n++) {
      Index i = 0;
      for (i = 0; i < struc_mol.size(); i++) {
        if (struc_mol[i] == m_mol_names[n]) {
          _add_rule(std::vector<Index>({i}));
          break;
        }
      }
      if (i == struc_mol.size())
        throw std::runtime_error(
            std::string("Format tag: '") + name() + "(" + m_mol_names[n] +
            ")' does not correspond to a viable composition.\n");
    }
  }
  return true;
}

/// \brief col_header returns: {'name(Au)', 'name(Pt)', ...}
std::vector<std::string> MolDependent::col_header(
    const Configuration &_tmplt) const {
  std::vector<std::string> col;
  for (Index c = 0; c < m_mol_names.size(); c++) {
    col.push_back(name() + "(" + m_mol_names[c] + ")");
  }
  return col;
}
}  // namespace ConfigIO_impl

namespace ConfigIO {

// --- Comp implementations -----------

const std::string Comp::Name = "comp";

const std::string Comp::Desc =
    "Parametric composition parameters, individual label as argument. "
    "Without argument, all values are printed. Ex: comp(a), comp(b), etc.";

/// \brief Returns the parametric composition
Eigen::VectorXd Comp::evaluate(const Configuration &config) const {
  return comp(config);
}

/// \brief Returns true if the PrimClex has composition axes
bool Comp::validate(const Configuration &config) const {
  return config.primclex().has_composition_axes();
}

/// \brief Expects arguments of the form 'comp(a)', 'comp(b)', etc.
bool Comp::parse_args(const std::string &args) {
  if (args.size() == 1) {
    _add_rule(std::vector<Index>({(Index)(args[0] - 'a')}));
  } else if (args.size() > 1) {
    throw std::runtime_error(std::string("Format tag: 'comp(") + args +
                             ") is invalid.\n");
    return false;
  }
  return true;
}

/// \brief col_header returns: {'comp(a)', 'comp(b)', ...}
std::vector<std::string> Comp::col_header(const Configuration &_tmplt) const {
  std::vector<std::string> col;
  for (Index c = 0; c < _index_rules().size(); c++) {
    col.push_back(name() + "(" + (char)('a' + _index_rules()[c][0]) + ")");
  }
  return col;
}

// --- CompN implementations -----------

const std::string CompN::Name = "comp_n";

const std::string CompN::Desc =
    "Number of each species per unit cell, including vacancies. "
    "No argument prints all available values. Ex: comp_n, comp_n(Au), "
    "comp_n(Pt), etc.";

/// \brief Returns the number of each species per unit cell
Eigen::VectorXd CompN::evaluate(const Configuration &config) const {
  return comp_n(config);
}

// --- SiteFrac implementations -----------

const std::string SiteFrac::Name = "site_frac";

const std::string SiteFrac::Desc =
    "Fraction of sites occupied by a species, including vacancies. "
    "No argument prints all available values. Ex: site_frac(Au), "
    "site_frac(Pt), etc.";

/// \brief Returns the site fraction
Eigen::VectorXd SiteFrac::evaluate(const Configuration &config) const {
  return site_frac(config);
}

// --- AtomFrac implementations -----------

const std::string AtomFrac::Name = "atom_frac";

const std::string AtomFrac::Desc =
    "Fraction of atoms that are a particular species, excluding vacancies.  "
    "Without argument, all values are printed. Ex: atom_frac(Au), "
    "atom_frac(Pt), etc.";

/// \brief Returns the site fraction
Eigen::VectorXd AtomFrac::evaluate(const Configuration &config) const {
  return species_frac(config);
}

// --- Corr implementations -----------

const std::string Corr::Name = "corr";

const std::string Corr::Desc =
    "Correlation values (evaluated basis functions, normalized per primitive "
    "cell). "
    "If no arguments, prints all correlations, using the basis set for the "
    "default "
    "cluster expansion as listed by 'casm settings -l'. "
    "If one argument, accepts either: "
    "1) a cluster expansion name, for example 'corr(formation_energy)', and "
    "evaluates all basis functions, or "
    "2) an integer index or range of indices of basis functions to evaluate, "
    "for example 'corr(6)', or 'corr(0:6)'. "
    "If two arguments, accepts cluster expansion name and an integer index or "
    "range of basis functions to evaluate, for example "
    "'corr(formation_energy,6)' "
    "or 'corr(formation_energy,0:6)'.";

/// \brief Returns the atom fraction
Eigen::VectorXd Corr::evaluate(const Configuration &config) const {
  Eigen::VectorXd corr;
  restricted_extensive_correlations(
      corr, config.configdof(), config.supercell().nlist(), m_clexulator,
      m_correlation_indices.data(), end_ptr(m_correlation_indices));
  corr /= ((double)config.supercell().volume());
  return corr;
}

/// \brief If not yet initialized, use the default clexulator from the PrimClex
bool Corr::init(const Configuration &_tmplt) const {
  if (!m_clexulator.initialized()) {
    const PrimClex &primclex = _tmplt.primclex();
    ClexDescription desc = m_clex_name.empty()
                               ? primclex.settings().default_clex()
                               : primclex.settings().clex(m_clex_name);
    m_clexulator = primclex.clexulator(desc.bset);
  }

  VectorXdAttribute<Configuration>::init(_tmplt);
  m_correlation_indices.clear();
  for (Index i = 0; i < _index_rules().size(); i++) {
    m_correlation_indices.push_back(_index_rules()[i][0]);
  }
  return true;
}

/// \brief Expects 'corr', 'corr(clex_name)', 'corr(index_expression)', or
/// 'corr(clex_name,index_expression)'
bool Corr::parse_args(const std::string &args) {
  std::vector<std::string> splt_vec;
  boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);

  if (!splt_vec.size()) {
    return true;
  } else if (splt_vec.size() == 1) {
    if ((splt_vec[0].find_first_not_of("0123456789") == std::string::npos) ||
        (splt_vec[0].find(':') != std::string::npos)) {
      _parse_index_expression(splt_vec[0]);
    } else {
      m_clex_name = splt_vec[0];
    }
  } else if (splt_vec.size() == 2) {
    m_clex_name = splt_vec[0];
    _parse_index_expression(splt_vec[1]);
  } else {
    std::stringstream ss;
    ss << "Too many arguments for 'corr'.  Received: " << args << "\n";
    throw std::runtime_error(ss.str());
  }
  return true;
}

// --- CorrContribution implementations -----------

const std::string CorrContribution::Name = "corr_contribution";

const std::string CorrContribution::Desc =
    "Correlation values (evaluated basis functions for a single unit cell, not "
    " normalized). The first argument is the linear unit cell index [0, "
    " scel_vol). The remaining arguments follow the same conventions as "
    "`corr` accepting any of `corr_contribution(linear_unitcell_index)` or "
    "`corr_contribution(linear_unitcell_index, clex_name)` or "
    "`corr_contribution(linear_unitcell_index, indices)`, or "
    "`corr_contribution(linear_unitcell_index, clex_name, indices)`. "
    "Coordinates of unit cells can be obtained from the `unitcells` "
    "information of the `info -m SupercellInfo` command.";

/// \brief Returns the atom fraction
Eigen::VectorXd CorrContribution::evaluate(const Configuration &config) const {
  return corr_contribution(m_linear_unitcell_index, config, m_clexulator);
}

/// \brief If not yet initialized, use the default clexulator from the PrimClex
bool CorrContribution::init(const Configuration &_tmplt) const {
  if (!m_clexulator.initialized()) {
    const PrimClex &primclex = _tmplt.primclex();
    ClexDescription desc = m_clex_name.empty()
                               ? primclex.settings().default_clex()
                               : primclex.settings().clex(m_clex_name);
    m_clexulator = primclex.clexulator(desc.bset);
  }

  VectorXdAttribute<Configuration>::init(_tmplt);
  return true;
}

///
/// Expects one of:
/// - 'corr_contribution(linear_unitcell_index)'
/// - 'corr(linear_unitcell_index), clex_name)'
/// - 'corr(linear_unitcell_index, index_expression)'
/// - 'corr(clex_name,index_expression)'
bool CorrContribution::parse_args(const std::string &args) {
  std::vector<std::string> splt_vec;
  boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);

  if (!splt_vec.size()) {
    return false;
  }
  if (splt_vec.size() == 1) {
    m_linear_unitcell_index = std::stol(splt_vec[0]);
  } else if (splt_vec.size() == 2) {
    m_linear_unitcell_index = std::stol(splt_vec[0]);
    if ((splt_vec[1].find_first_not_of("0123456789") == std::string::npos) ||
        (splt_vec[1].find(':') != std::string::npos)) {
      _parse_index_expression(splt_vec[1]);
    } else {
      m_clex_name = splt_vec[1];
    }
  } else if (splt_vec.size() == 3) {
    m_linear_unitcell_index = std::stol(splt_vec[0]);
    m_clex_name = splt_vec[1];
    _parse_index_expression(splt_vec[2]);
  } else {
    std::stringstream ss;
    ss << "Too many arguments for 'corr_contribution'.  Received: " << args
       << "\n";
    throw std::runtime_error(ss.str());
  }
  return true;
}

// --- PointCorr implementations -----------

const std::string PointCorr::Name = "point_corr";

const std::string PointCorr::Desc =
    "Point correlation values (evaluated basis functions for a single site, "
    "normalized by cluster orbit size). The first argument is the linear unit "
    "cell index [0, scel_vol). The second argument is the index of the site in "
    "the neighbor list. For periodic cluster functions, the `neighbor_index` "
    "counts over sites in the primitive cell that have site degrees of "
    "freedom. If all sites in the prim have site degrees of freedom, then the "
    "`neighbor_index` is equal to the `sublattice_index` of that site in the "
    "prim. For local cluster functions, the `neighbor_index` counts over sites "
    "in the neighborhood that have site degrees of freedom. The remaining "
    "arguments follow the same conventions as `corr` accepting any of "
    "`point_corr(linear_unitcell_index, neighbor_index)` or "
    "`point_corr(linear_unitcell_index, neighbor_index, clex_name)` or "
    "`point_corr(linear_unitcell_index, neighbor_index, indices)`, or "
    "`point_corr(linear_unitcell_index, neighbor_index, clex_name, "
    "indices)`. Coordinates of sites can be obtained from the `info -m "
    "NeighborListInfo` method.";

Eigen::VectorXd PointCorr::evaluate(const Configuration &config) const {
  return point_corr(m_linear_unitcell_index, m_neighbor_index, config,
                    m_clexulator);
}

/// \brief If not yet initialized, use the default clexulator from the PrimClex
bool PointCorr::init(const Configuration &_tmplt) const {
  if (!m_clexulator.initialized()) {
    const PrimClex &primclex = _tmplt.primclex();
    ClexDescription desc = m_clex_name.empty()
                               ? primclex.settings().default_clex()
                               : primclex.settings().clex(m_clex_name);
    m_clexulator = primclex.clexulator(desc.bset);
  }

  VectorXdAttribute<Configuration>::init(_tmplt);
  return true;
}

///
/// Expects one of:
/// - 'point_corr(linear_unitcell_index, neighbor_index)'
/// - 'point_corr(linear_unitcell_index, neighbor_index, clex_name)'
/// - 'point_corr(linear_unitcell_index, neighbor_index, index_expression)'
/// - 'point_corr(linear_unitcell_index, neighbor_index, clex_name,
/// index_expression)'
bool PointCorr::parse_args(const std::string &args) {
  std::vector<std::string> splt_vec;
  boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);

  if (splt_vec.size() < 2) {
    return false;
  }
  if (splt_vec.size() == 2) {
    m_linear_unitcell_index = std::stol(splt_vec[0]);
    m_neighbor_index = std::stol(splt_vec[1]);
  } else if (splt_vec.size() == 3) {
    m_linear_unitcell_index = std::stol(splt_vec[0]);
    m_neighbor_index = std::stol(splt_vec[1]);
    if ((splt_vec[2].find_first_not_of("0123456789") == std::string::npos) ||
        (splt_vec[2].find(':') != std::string::npos)) {
      _parse_index_expression(splt_vec[2]);
    } else {
      m_clex_name = splt_vec[2];
    }
  } else if (splt_vec.size() == 4) {
    m_linear_unitcell_index = std::stol(splt_vec[0]);
    m_neighbor_index = std::stol(splt_vec[1]);
    m_clex_name = splt_vec[2];
    _parse_index_expression(splt_vec[3]);
  } else {
    std::stringstream ss;
    ss << "Too many arguments for 'point_corr'.  Received: " << args << "\n";
    throw std::runtime_error(ss.str());
  }
  return true;
}

// --- AllCorrContribution implementations -----------

const std::string AllCorrContribution::Name = "all_corr_contribution";

const std::string AllCorrContribution::Desc =
    "Correlation values (evaluated basis functions for a single unit cell, not "
    "normalized), for every unitcell in the supercell. The output is a matrix, "
    "with each row corresponding to a unitcell. The mean over rows is equal to "
    "the `corr` output. The arguments follow the same conventions as `corr` "
    "accepting any of:\n\n"
    "- `all_corr_contribution`\n\n"
    "- `all_corr_contribution(clex_name)`\n\n"
    "- `all_corr_contribution(indices)`\n\n"
    "- `all_corr_contribution(clex_name, indices)`. \n\n"
    "Coordinates of unit cells can be obtained from the `unitcells` "
    "information of the `info -m SupercellInfo` command.";

Eigen::MatrixXd AllCorrContribution::evaluate(
    const Configuration &config) const {
  return all_corr_contribution(config, m_clexulator);
}

/// \brief If not yet initialized, use the default clexulator from the PrimClex
bool AllCorrContribution::init(const Configuration &_tmplt) const {
  if (!m_clexulator.initialized()) {
    const PrimClex &primclex = _tmplt.primclex();
    ClexDescription desc = m_clex_name.empty()
                               ? primclex.settings().default_clex()
                               : primclex.settings().clex(m_clex_name);
    m_clexulator = primclex.clexulator(desc.bset);
  }

  MatrixXdAttribute<Configuration>::init(_tmplt);
  return true;
}

///
/// Expects one of:
/// - 'all_corr_contribution'
/// - 'all_corr_contribution(clex_name)'
/// - 'all_corr_contribution(index_expression)'
/// - 'all_corr_contribution(clex_name, index_expression)'
bool AllCorrContribution::parse_args(const std::string &args) {
  std::vector<std::string> splt_vec;
  boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);

  if (!splt_vec.size()) {
    return true;
  }
  if (splt_vec.size() == 1) {
    if ((splt_vec[0].find_first_not_of("0123456789") == std::string::npos) ||
        (splt_vec[0].find(':') != std::string::npos)) {
      _parse_index_expression(splt_vec[0]);
    } else {
      m_clex_name = splt_vec[0];
    }
  } else if (splt_vec.size() == 2) {
    m_clex_name = splt_vec[0];
    _parse_index_expression(splt_vec[1]);
  } else {
    std::stringstream ss;
    ss << "Too many arguments for 'all_corr_contribution'.  Received: " << args
       << "\n";
    throw std::runtime_error(ss.str());
  }
  return true;
}

// --- SiteCentricCorrelations implementations -----------

const std::string SiteCentricCorrelations::Name = "site_centric_correlations";

const std::string SiteCentricCorrelations::Desc =
    "Site-centric correlation values, the sum of all cluster functions that "
    "include a particular site, normalized by cluster orbit size. The output "
    "is a JSON object, with a \"value\" matrix of size (n_sites, "
    "n_correlations) and an \"asymmetric_unit_indices\" array of size n_sites. "
    "Each row in the \"value\" matrix contains the point correlations for a "
    "site. The sum over rows normalized by supercell size (number of unit "
    "cells) is equal to the global mean correlations (`corr`) times the "
    "cluster size for functions associated with a given column (due to "
    "multiple counting of cluster functions for each site involved in the "
    "cluster). The i-th element of \"asymmetric_unit_indices\" gives the "
    "asymmetric unit index for the site's sublattice (equal indices indicates "
    "sites on symmetrically equivalent sublattices). Accepts either zero "
    "arguments (`site_centric_correlations`), and uses the default cluster "
    "expansion or one argument, (`site_centric_correlations(clex_name)`) and "
    "uses the named cluster expansion. More detailed coordinates of sites for "
    "any supercell can be obtained from the `info -m NeighborListInfo` method.";

jsonParser SiteCentricCorrelations::evaluate(
    const Configuration &config) const {
  Supercell const &supercell = config.supercell();
  auto const &origin_unitcellcoords =
      m_neighborhood_info->point_corr_unitcellcoord;

  std::vector<Index> asymmetric_unit_indices;
  for (xtal::UnitCellCoord unitcellcoord : origin_unitcellcoords) {
    for (Index v = 0; v < supercell.volume(); ++v) {
      asymmetric_unit_indices.push_back(
          m_sublattice_to_asymmetric_unit[unitcellcoord.sublattice()]);
    }
  }

  jsonParser json;
  json["asymmetric_unit_indices"] = asymmetric_unit_indices;
  json["value"] = all_point_corr(config, m_clexulator);
  return json;
}

bool SiteCentricCorrelations::validate(const Configuration &config) const {
  if (!m_neighborhood_info) {
    CASM::err_log() << "Error in \"site_centric_correlations\": not properly "
                       "initialized (unknown error)."
                    << std::endl;
    return false;
  }
  return true;
}

/// \brief If not yet initialized, use the default clexulator from the PrimClex
bool SiteCentricCorrelations::init(const Configuration &_tmplt) const {
  if (!m_clexulator.initialized()) {
    const PrimClex &primclex = _tmplt.primclex();
    auto const &prim = primclex.prim();
    ClexDescription desc = m_clex_name.empty()
                               ? primclex.settings().default_clex()
                               : primclex.settings().clex(m_clex_name);
    ClexBasisSpecs const &basis_set_specs = primclex.basis_set_specs(desc.bset);
    ClusterSpecs const &cluster_specs = *basis_set_specs.cluster_specs;
    if (cluster_specs.periodicity_type() !=
        CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC) {
      throw std::runtime_error(
          "Error: \"site_centric_correlations\" is only valid for periodic "
          "cluster expansions.");
    }
    if (static_cast<PeriodicMaxLengthClusterSpecs const &>(cluster_specs)
            .generating_group.size() != prim.factor_group().size()) {
      throw std::runtime_error(
          "Error: \"site_centric_correlations\" is only valid when the "
          "generating group is the prim factor group.");
    }

    m_clexulator = primclex.clexulator(desc.bset);
    m_neighborhood_info = &primclex.neighborhood_info(desc.bset);

    adapter::Adapter<xtal::SymOpVector, SymGroup> adapter;
    // set of [sets of basis indices of equivalent sites]
    std::set<std::set<Index>> orbits =
        xtal::make_asymmetric_unit(prim, adapter(prim.factor_group()));

    m_sublattice_to_asymmetric_unit.resize(prim.basis().size());
    Index asym_unit_index = 0;
    for (std::set<Index> const &orbit : orbits) {
      for (Index const &sublat : orbit) {
        m_sublattice_to_asymmetric_unit[sublat] = asym_unit_index;
      }
      asym_unit_index++;
    }
  }

  BaseValueFormatter<jsonParser, Configuration>::init(_tmplt);
  return true;
}

/// Expects one of:
/// - 'site_centric_correlations'
/// - 'site_centric_correlations(clex_name)'
bool SiteCentricCorrelations::parse_args(const std::string &args) {
  std::vector<std::string> splt_vec;
  boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);

  if (!splt_vec.size()) {
    return true;
  } else if (splt_vec.size() == 1) {
    m_clex_name = splt_vec[0];
  } else {
    std::stringstream ss;
    ss << "Too many arguments for 'site_centric_correlations'.  Received: "
       << args << "\n";
    throw std::runtime_error(ss.str());
  }
  return true;
}

// --- AllPointCorr implementations -----------

const std::string AllPointCorr::Name = "all_point_corr";

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

const std::string AllPointCorr::Desc =
    "Point correlation values, the sum of all cluster functions that include a "
    "particular site, normalized by cluster orbit size, output as a JSON "
    "object. The \"value\" is a matrix of size (n_sites, n_correlations). Each "
    "row in the \"value\" matrix contains the point correlations for one site. "
    "The interpretation and application differs for periodic vs local cluster "
    "expansions. The sum over rows normalized by supercell size (number of "
    "unit cells) is equal to the global mean correlations (`corr`) times the "
    "cluster size for functions associated with a given column (due to "
    "multiple counting of cluster functions for each site involved in the "
    "cluster). Information about the sites is provided by the arrays (each of "
    "size n_sites) \"asymmetric_unit_indices\" (equivalent indices indicates "
    "symmetrically equivalent sites), \"linear_unitcell_index\" (indicates "
    "which unit cell the functions were evaluated for), \"unitcellcoord\" (the "
    "integral site coordinates (b, i, j, k) of the site, where b is the "
    "sublattice index and (i,j,k) are integral coordinates of the unit cell "
    "containing the site). Accepts either zero arguments (`all_point_corr`), "
    "and uses the default cluster expansion or one argument, "
    "(`all_point_corr(clex_name)`) and uses the named cluster expansion.";

jsonParser AllPointCorr::evaluate(const Configuration &config) const {
  Supercell const &supercell = config.supercell();
  jsonParser json;
  json["asymmetric_unit_indices"] = make_all_point_corr_asymmetric_unit_indices(
      *m_neighborhood_info, supercell.sym_info());
  json["linear_unitcell_indices"] = make_all_point_corr_linear_unitcell_indices(
      *m_neighborhood_info, supercell.sym_info());
  json["unitcellcoord"] = make_all_point_corr_unitcellcoord(
      *m_neighborhood_info, supercell.sym_info());
  json["value"] = all_point_corr(config, m_clexulator);
  return json;
}

bool AllPointCorr::validate(const Configuration &config) const {
  if (!m_neighborhood_info) {
    CASM::err_log() << "Error in AllPointCorr: not properly initialized."
                    << std::endl;
    return false;
  }
  return true;
}

/// \brief If not yet initialized, use the default clexulator from the PrimClex
bool AllPointCorr::init(const Configuration &_tmplt) const {
  if (!m_clexulator.initialized()) {
    const PrimClex &primclex = _tmplt.primclex();
    ClexDescription desc = m_clex_name.empty()
                               ? primclex.settings().default_clex()
                               : primclex.settings().clex(m_clex_name);
    m_clexulator = primclex.clexulator(desc.bset);
    m_neighborhood_info = &primclex.neighborhood_info(desc.bset);
  }

  BaseValueFormatter<jsonParser, Configuration>::init(_tmplt);
  return true;
}

///
/// Expects one of:
/// - 'all_point_corr'
/// - 'all_point_corr(clex_name)'
bool AllPointCorr::parse_args(const std::string &args) {
  std::vector<std::string> splt_vec;
  boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);

  if (!splt_vec.size()) {
    return true;
  } else if (splt_vec.size() == 1) {
    m_clex_name = splt_vec[0];
  } else {
    std::stringstream ss;
    ss << "Too many arguments for 'all_point_corr'.  Received: " << args
       << "\n";
    throw std::runtime_error(ss.str());
  }
  return true;
}

// --- GradCorr implementations -----------

const std::string GradCorr::Name = "gradcorr";

const std::string GradCorr::Desc =
    "Gradiant of correlation values (evaluated basis functions), with respect "
    "to a specified "
    "degree of freedom (DoF). For each configuration, output is a (D*N x C) "
    "matrix, where 'D' "
    "is DoF dimension, 'N' is either 1 (for global DoF) or number of sites in "
    "the configuration "
    "(for site DoF), and 'C' is number of basis functions. Gradient components "
    "are ordered such "
    "that components corresponding to particular site are listed in "
    "consecutive rows. Requires "
    "at least one argument, specifying the DoF with repect to which gradient "
    "is taken [e.g., 'gradcorr(disp)']. "
    "Basis functions are the basis set for the default cluster expansion, as "
    "listed by 'casm settings -l', "
    "unless otherwise specified. Accepts up to three additional arguments: "
    "1) a cluster expansion name, e.g. 'gradcorr(disp,formation_energy)' "
    "and/or "
    "2) a pair of indices, or index ranges, e.g. 'gradcorr(disp,5,2)', "
    "'gradcorr(disp,:,3:5)', 'gradcorr(disp,0:4,3)'";

/// \brief Returns the atom fraction
Eigen::MatrixXd GradCorr::evaluate(const Configuration &config) const {
  return gradcorrelations(config, m_clexulator, m_key);
}

/// \brief If not yet initialized, use the default clexulator from the PrimClex
bool GradCorr::init(const Configuration &_tmplt) const {
  if (!m_clexulator.initialized()) {
    const PrimClex &primclex = _tmplt.primclex();
    ClexDescription desc = m_clex_name.empty()
                               ? primclex.settings().default_clex()
                               : primclex.settings().clex(m_clex_name);
    m_clexulator = primclex.clexulator(desc.bset);
  }

  MatrixXdAttribute<Configuration>::init(_tmplt);
  return true;
}

/// \brief Expects 'corr', 'corr(clex_name)', 'corr(index_expression)', or
/// 'corr(clex_name,index_expression)'
bool GradCorr::parse_args(const std::string &args) {
  // std::cout << "parsing args: " << args << "\n";
  std::vector<std::string> split_vec;
  boost::split(split_vec, args, boost::is_any_of(","),
               boost::token_compress_on);
  // std::cout << "after split: " << split_vec << "\n";
  if (!split_vec.size()) {
    throw std::runtime_error(
        "'gradcorr' query requires at least one argument, corresponding to the "
        "independent variable wrt which gradient is to be computed.");
    return false;
  } else if (split_vec.size() > 4) {
    std::stringstream ss;
    ss << "Too many arguments for 'gradcorr'.  Received: " << args << "\n";
    throw std::runtime_error(ss.str());
  }

  boost::erase_all(split_vec[0], "'");
  // std::cout << "Now split_vec[0] is " << split_vec[0] << "\n";
  m_key = split_vec[0];

  for (Index i = 1; i < split_vec.size(); ++i) {
    if ((split_vec[i].find_first_not_of("0123456789") == std::string::npos) ||
        (split_vec[i].find(':') != std::string::npos)) {
      _parse_index_expression(split_vec[i] + "," + split_vec[i + 1]);
      ++i;
    } else {
      m_clex_name = split_vec[i];
    }
  }
  return true;
}

// --- Clex implementations -----------

const std::string Clex::Name = "clex";

const std::string Clex::Desc =
    "Predicted property value."
    " Accepts arguments ($clex_name,$norm)."
    " ($clex_name is a cluster expansion name as listed by 'casm settings -l', "
    "default=the default clex)"
    " ($norm is the normalization, either 'per_species', or 'per_unitcell' "
    "<--default)";

Clex::Clex() : ScalarAttribute<Configuration>(Name, Desc) { parse_args(""); }

Clex::Clex(const Clexulator &clexulator, const ECIContainer &eci,
           const Norm<Configuration> &norm)
    : ScalarAttribute<Configuration>(Name, Desc),
      m_clexulator(clexulator),
      m_eci(eci),
      m_norm(norm.clone()) {}

/// \brief Returns the atom fraction
double Clex::evaluate(const Configuration &config) const {
  return m_eci * correlations(config, m_clexulator) / _norm(config);
}

/// \brief Clone using copy constructor
std::unique_ptr<Clex> Clex::clone() const {
  return std::unique_ptr<Clex>(this->_clone());
}

/// \brief If not yet initialized, use the default cluster expansion from the
/// PrimClex
bool Clex::init(const Configuration &_tmplt) const {
  if (!m_clexulator.initialized()) {
    const PrimClex &primclex = _tmplt.primclex();
    ClexDescription desc = m_clex_name.empty()
                               ? primclex.settings().default_clex()
                               : primclex.settings().clex(m_clex_name);
    m_clexulator = primclex.clexulator(desc.bset);
    m_eci = primclex.eci(desc);
    if (m_eci.index().back() >= m_clexulator.corr_size()) {
      Log &err_log = CASM::err_log();
      err_log.error<Log::standard>("bset and eci mismatch");
      err_log << "basis set size: " << m_clexulator.corr_size() << std::endl;
      err_log << "max eci index: " << m_eci.index().back() << std::endl;
      throw std::runtime_error("Error: bset and eci mismatch");
    }
  }
  return true;
}

/// \brief Expects 'clex', 'clex(formation_energy)', or
/// 'clex(formation_energy,per_species)'
bool Clex::parse_args(const std::string &args) {
  std::vector<std::string> splt_vec;
  boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);

  m_clex_name = "";
  if (splt_vec.size()) {
    m_clex_name = splt_vec[0];
  }

  m_norm = notstd::make_cloneable<Norm<Configuration>>();
  if (splt_vec.size() == 2) {
    if (splt_vec[1] == "per_unitcell") {
      m_norm = notstd::make_cloneable<Norm<Configuration>>();
    } else if (splt_vec[1] == "per_species") {
      m_norm = notstd::make_cloneable<NormPerSpecies>();
    } else {
      std::stringstream ss;
      ss << "Error parsing second argument for 'clex'.  Received: " << args
         << "\n";
      throw std::runtime_error(ss.str());
    }
  }

  if (splt_vec.size() > 2) {
    std::stringstream ss;
    ss << "Too many arguments for 'clex'.  Received: " << args << "\n";
    throw std::runtime_error(ss.str());
  }

  return true;
}

/// \brief Returns the normalization
double Clex::_norm(const Configuration &config) const {
  return (*m_norm)(config);
}

/// \brief Clone using copy constructor
Clex *Clex::_clone() const { return new Clex(*this); }

/*End ConfigIO*/
}  // namespace ConfigIO

namespace ConfigIO {

GenericConfigFormatter<std::string> configname() {
  return GenericConfigFormatter<std::string>(
      "configname", "Configuration name, in the form 'SCEL#_#_#_#_#_#_#/#'",
      [](const Configuration &config) -> std::string { return config.name(); });
}

GenericConfigFormatter<std::string> scelname() {
  return GenericConfigFormatter<std::string>(
      "scelname", "Supercell name, in the form 'SCEL#_#_#_#_#_#_#'",
      [](const Configuration &config) -> std::string {
        return config.supercell().name();
      });
}

GenericConfigFormatter<std::string> calc_status() {
  return GenericConfigFormatter<std::string>(
      "calc_status",
      "Status of calculation."
      "A calculation that has been run successfully will be marked 'complete'."
      "Mapped or imported configurations may have 'is_calculated = 1' without "
      "'calc_status = complete'.",
      [](const Configuration &config) -> std::string {
        return CASM::calc_status<Configuration>(config);
      },
      [](const Configuration &config) -> bool {
        return CASM::has_calc_status<Configuration>(config);
      });
}

GenericConfigFormatter<std::string> failure_type() {
  return GenericConfigFormatter<std::string>(
      "failure_type", "Reason for calculation failure.",
      [](const Configuration &config) -> std::string {
        return CASM::failure_type<Configuration>(config);
      },
      [](const Configuration &config) -> bool {
        return CASM::has_failure_type<Configuration>(config);
      });
}

GenericConfigFormatter<Index> scel_size() {
  return GenericConfigFormatter<Index>(
      "scel_size",
      "Supercell volume, given as the integer number of primitive cells",
      [](const Configuration &config) -> Index {
        return config.supercell().volume();
      });
}

GenericConfigFormatter<Index> multiplicity() {
  return GenericConfigFormatter<Index>(
      "multiplicity",
      "Symmetric multiplicity of the configuration, excluding translational "
      "equivalents.",
      [](const Configuration &config) -> Index {
        return config.multiplicity();
      });
}

GenericConfigFormatter<std::string> point_group_name() {
  return GenericConfigFormatter<std::string>(
      "point_group_name", "Name of the configuration's point group.",
      [](const Configuration &config) -> std::string {
        return config.point_group_name();
      });
}

// deprecated for 'energy'
GenericConfigFormatter<double> relaxed_energy() {
  return GenericConfigFormatter<double>(
      "relaxed_energy",
      "DFT energy, normalized per primitive cell (deprecated for `energy`)",
      CASM::energy, has_energy);
}

GenericConfigFormatter<double> energy() {
  return GenericConfigFormatter<double>(
      "energy", "DFT energy, normalized per primitive cell", CASM::energy,
      has_energy);
}

// deprecated for 'energy_per_species'
GenericConfigFormatter<double> relaxed_energy_per_species() {
  return GenericConfigFormatter<double>(
      "relaxed_energy_per_atom",
      "DFT energy, normalized per atom (deprecated for `energy_per_atom`)",
      CASM::energy_per_species, has_energy);
}

GenericConfigFormatter<double> energy_per_species() {
  return GenericConfigFormatter<double>("energy_per_atom",
                                        "DFT energy, normalized per atom",
                                        CASM::energy_per_species, has_energy);
}

GenericConfigFormatter<double> reference_energy() {
  return GenericConfigFormatter<double>(
      "reference_energy",
      "reference energy, normalized per primitive cell, as determined by "
      "current reference states",
      CASM::reference_energy, has_reference_energy);
}

GenericConfigFormatter<double> reference_energy_per_species() {
  return GenericConfigFormatter<double>(
      "reference_energy_per_atom",
      "reference energy, normalized per atom, as determined by current "
      "reference states",
      CASM::reference_energy_per_species, has_reference_energy);
}

GenericConfigFormatter<double> formation_energy() {
  return GenericConfigFormatter<double>(
      "formation_energy",
      "DFT formation energy, normalized per primitive cell and measured "
      "relative to current reference states",
      CASM::formation_energy, has_formation_energy);
}

GenericConfigFormatter<double> formation_energy_per_species() {
  return GenericConfigFormatter<double>(
      "formation_energy_per_atom",
      "DFT formation energy, normalized per atom and measured relative to "
      "current reference states",
      CASM::formation_energy_per_species, has_formation_energy);
}

/*Generic1DDatumFormatter<std::vector<double>, Configuration
  >relaxation_strain() { return Generic1DDatumFormatter<std::vector<double>,
  Configuration >("relaxation_strain", "Green-Lagrange strain of dft-relaxed
  configuration, relative to the ideal crystal.  Ordered as [E(0,0), E(1,1),
  E(2,2), E(1,2), E(0,2), E(0,1)].  Accepts index as argument on interval
  [0,5]", CASM::relaxation_strain, has_relaxation_strain,
  [](const std::vector<double> &cont)->Index{
    return 6;
  });
  }*/

GenericConfigFormatter<bool> is_calculated() {
  return GenericConfigFormatter<bool>(
      "is_calculated",
      "True (1) if all current properties have been been calculated for the "
      "configuration",
      [](const Configuration &config) -> bool {
        return CASM::is_calculated(config);
      });
}

GenericConfigFormatter<bool> is_primitive() {
  return GenericConfigFormatter<bool>("is_primitive",
                                      "True (1) if the configuration cannot be "
                                      "described within a smaller supercell",
                                      CASM::is_primitive);
}

GenericConfigFormatter<bool> is_canonical() {
  return GenericConfigFormatter<bool>(
      "is_canonical",
      "True (1) if the configuration cannot be transformed by symmetry to a "
      "configuration with higher lexicographic order",
      CASM::is_canonical);
}

GenericConfigFormatter<double> rms_force() {
  return GenericConfigFormatter<double>(
      "rms_force",
      "Root-mean-square forces of relaxed configurations, determined from DFT "
      "(eV/Angstr.)",
      CASM::rms_force, has_rms_force);
}

GenericConfigFormatter<double> atomic_deformation() {
  return GenericConfigFormatter<double>(
      "atomic_deformation",
      "Cost function that describes the degree to which basis sites have "
      "relaxed",
      CASM::atomic_deformation, has_atomic_deformation);
}

GenericConfigFormatter<double> lattice_deformation() {
  return GenericConfigFormatter<double>(
      "lattice_deformation",
      "Cost function that describes the degree to which lattice has relaxed.",
      CASM::lattice_deformation, has_lattice_deformation);
}

GenericConfigFormatter<double> volume_relaxation() {
  return GenericConfigFormatter<double>(
      "volume_relaxation",
      "Change in volume due to relaxation, expressed as the ratio V/V_0.",
      CASM::volume_relaxation, has_volume_relaxation);
}

GenericConfigFormatter<double> relaxed_magmom() {
  return GenericConfigFormatter<double>(
      "relaxed_magmom",
      "Relaxed magnetic moment, normalized per primative cell.",
      CASM::relaxed_magmom, has_relaxed_magmom);
}

GenericConfigFormatter<double> relaxed_magmom_per_species() {
  return GenericConfigFormatter<double>(
      "relaxed_magmom_per_atom",
      "Relaxed magnetic moment, normalized per atom.",
      CASM::relaxed_magmom_per_species, has_relaxed_magmom);
}

ConfigIO::GenericConfigFormatter<jsonParser> structure() {
  return GenericConfigFormatter<jsonParser>(
      "structure",
      "Structure resulting from application of DoF, formatted as JSON",
      [](Configuration const &configuration) {
        jsonParser json = jsonParser::object();
        to_json(make_simple_structure(configuration), json);
        return json;
      });
}

ConfigIO::GenericConfigFormatter<jsonParser> structure_with_vacancies() {
  return GenericConfigFormatter<jsonParser>(
      "structure_with_vacancies",
      "Structure resulting from application of DoF, including vacancies, "
      "formatted as JSON",
      [](Configuration const &configuration) {
        jsonParser json = jsonParser::object();
        std::set<std::string> const &excluded_species = {};
        to_json(make_simple_structure(configuration), json, excluded_species);
        return json;
      });
}

ConfigIO::GenericConfigFormatter<jsonParser> config() {
  return GenericConfigFormatter<jsonParser>(
      "config", "All degrees of freedom (DoF), formatted as JSON",
      [](Configuration const &configuration) {
        jsonParser json = jsonParser::object();
        to_json(configuration, json);
        return json;
      });
}

ConfigIO::GenericConfigFormatter<jsonParser> properties() {
  return GenericConfigFormatter<jsonParser>(
      "properties",
      "Configuration properties, formatted as JSON. In case of multiple "
      "structures mapping to the same configuration, the conflict score is "
      "used to determine which properties are returned.",
      [](Configuration const &configuration) {
        return jsonParser{configuration.calc_properties()};
      });
}

ConfigIO::GenericConfigFormatter<jsonParser> all_mapped_properties() {
  return GenericConfigFormatter<jsonParser>(
      "all_mapped_properties",
      "All properties mapped to a configuration, for the current calctype, "
      "sorted by conflict score, and formatted as JSON.",
      [](Configuration const &configuration) {
        PrimClex const &primclex = configuration.primclex();
        std::string calctype = primclex.settings().default_clex().calctype;
        auto const &db = primclex.const_db_props<Configuration>(calctype);
        auto all_origins = db.all_origins(configuration.name());
        jsonParser json;
        json["conflict_score_method"] = db.score_method(configuration.name());
        json["mapped_properties"].put_array();
        for (std::string origin_name : all_origins) {
          jsonParser tjson;
          auto it = db.find_via_origin(origin_name);
          if (it != db.end()) {
            tjson["properties"] = *it;
            double conflict_score = db.score(origin_name);
            if (conflict_score != std::numeric_limits<double>::max()) {
              tjson["conflict_score"] = db.score(origin_name);
            } else {
              tjson["conflict_score"] = jsonParser::null();
            }
          } else {
            tjson["properties"] = jsonParser::null();
            tjson["conflict_score"] = jsonParser::null();
          }
          json["mapped_properties"].push_back(tjson);
        }
        return json;
      });
}

ConfigIO::GenericConfigFormatter<jsonParser> mapped_structure() {
  return GenericConfigFormatter<jsonParser>(
      "mapped_structure",
      "Mapped structure, formed by application of DoF and properties, and "
      "formatted as JSON. In case of multiple structures mapping to the same "
      "configuration, the conflict score is used to determine which properties "
      "are returned. The structure is written in the `mapped` orientation, a "
      "rigid transformation from the original`unmapped` orientation.",
      [](Configuration const &configuration) {
        jsonParser json = jsonParser::object();
        bool apply_properties = true;
        to_json(make_simple_structure(configuration, {}, apply_properties),
                json);
        return json;
      });
}

ConfigIO::GenericConfigFormatter<jsonParser> all_mapped_structures() {
  return GenericConfigFormatter<jsonParser>(
      "all_mapped_structures",
      "All structures that have been mapped to a configuration, formed by "
      "application of DoF and properties, sorted by conflict score, and "
      "formatted as JSON. The structures are written in the `mapped` "
      "orientation, a rigid transformation from their original `unmapped` "
      "orientations.",
      [](Configuration const &configuration) {
        PrimClex const &primclex = configuration.primclex();
        std::string calctype = primclex.settings().default_clex().calctype;
        auto const &db = primclex.const_db_props<Configuration>(calctype);
        auto all_origins = db.all_origins(configuration.name());
        jsonParser json;
        json["conflict_score_method"] = db.score_method(configuration.name());
        json["mapped_structures"].put_array();
        for (std::string origin_name : all_origins) {
          jsonParser tjson;
          tjson["origin"] = origin_name;
          auto it = db.find_via_origin(origin_name);
          if (it != db.end()) {
            jsonParser json = jsonParser::object();
            auto const &supercell = configuration.supercell();
            auto const &configdof = configuration.configdof();
            bool apply_properties = true;
            to_json(make_simple_structure(supercell, configdof, *it, {},
                                          apply_properties),
                    tjson["structure"]);
            double conflict_score = db.score(origin_name);
            if (conflict_score != std::numeric_limits<double>::max()) {
              tjson["conflict_score"] = db.score(origin_name);
            } else {
              tjson["conflict_score"] = jsonParser::null();
            }
          } else {
            tjson["structure"] = jsonParser::null();
            tjson["conflict_score"] = jsonParser::null();
          }
          json["mapped_structures"].push_back(tjson);
        }
        return json;
      });
}

ConfigIO::GenericConfigFormatter<std::string> poscar() {
  return GenericConfigFormatter<std::string>(
      "poscar",
      "Structure resulting from application of DoF, formatted as VASP POSCAR",
      [](Configuration const &configuration) {
        std::stringstream ss;
        std::string name;
        if (configuration.id() == "none") {
          name = configuration.supercell().name() + "/none";
        } else {
          name = configuration.name();
        }
        VaspIO::PrintPOSCAR p{make_simple_structure(configuration), name};
        p.sort();
        p.print(ss);
        return ss.str();
      });
}

ConfigIO::GenericConfigFormatter<std::string> poscar_with_vacancies() {
  return GenericConfigFormatter<std::string>(
      "poscar_with_vacancies",
      "Structure resulting from application of DoF, including vacancies, "
      "formatted as VASP POSCAR",
      [](Configuration const &configuration) {
        std::stringstream ss;
        std::string name;
        if (configuration.id() == "none") {
          name = configuration.supercell().name() + "/none";
        } else {
          name = configuration.name();
        }
        VaspIO::PrintPOSCAR p{make_simple_structure(configuration), name};
        p.ignore() = {};
        p.sort();
        p.print(ss);
        return ss.str();
      });
}

/*End ConfigIO*/
}  // namespace ConfigIO

template <>
StringAttributeDictionary<Configuration>
make_string_dictionary<Configuration>() {
  using namespace ConfigIO;
  StringAttributeDictionary<Configuration> dict;

  dict.insert(name<Configuration>(), configname(), alias<Configuration>(),
              alias_or_name<Configuration>(), scelname(), calc_status(),
              failure_type(), point_group_name(), poscar(),
              poscar_with_vacancies());

  return dict;
}

template <>
BooleanAttributeDictionary<Configuration>
make_boolean_dictionary<Configuration>() {
  using namespace ConfigIO;
  BooleanAttributeDictionary<Configuration> dict;

  dict.insert(is_calculated(), is_canonical(), is_primitive(),
              DB::Selected<Configuration>(), OnClexHull(), OnHull());

  return dict;
}

template <>
IntegerAttributeDictionary<Configuration>
make_integer_dictionary<Configuration>() {
  using namespace ConfigIO;
  IntegerAttributeDictionary<Configuration> dict;

  dict.insert(scel_size(), multiplicity());

  return dict;
}

template <>
ScalarAttributeDictionary<Configuration>
make_scalar_dictionary<Configuration>() {
  using namespace ConfigIO;
  ScalarAttributeDictionary<Configuration> dict;

  dict.insert(Clex(), HullDist(), ClexHullDist(), Novelty(), relaxed_energy(),
              energy(), relaxed_energy_per_species(), energy_per_species(),
              reference_energy(), reference_energy_per_species(),
              formation_energy(), formation_energy_per_species(), rms_force(),
              atomic_deformation(), lattice_deformation(), volume_relaxation(),
              relaxed_magmom(), relaxed_magmom_per_species());

  return dict;
}

template <>
VectorXiAttributeDictionary<Configuration>
make_vectorxi_dictionary<Configuration>() {
  using namespace ConfigIO;
  VectorXiAttributeDictionary<Configuration> dict;
  return dict;
}

template <>
VectorXdAttributeDictionary<Configuration>
make_vectorxd_dictionary<Configuration>() {
  using namespace ConfigIO;
  VectorXdAttributeDictionary<Configuration> dict;

  dict.insert(AtomFrac(), Comp(), CompN(), Corr(), CorrContribution(),
              PointCorr(), RelaxationStrain(), DoFStrain(), SiteFrac(),
              StrucScore());

  return dict;
}

template <>
MatrixXdAttributeDictionary<Configuration>
make_matrixxd_dictionary<Configuration>() {
  using namespace ConfigIO;
  MatrixXdAttributeDictionary<Configuration> dict;

  dict.insert(AllCorrContribution(), GradCorr());

  return dict;
}

template <>
DataFormatterDictionary<Configuration,
                        BaseValueFormatter<jsonParser, Configuration>>
make_json_dictionary<Configuration>() {
  using namespace ConfigIO;
  DataFormatterDictionary<Configuration,
                          BaseValueFormatter<jsonParser, Configuration>>
      dict;

  dict.insert(structure(), structure_with_vacancies(), config(), properties(),
              all_mapped_properties(), mapped_structure(),
              all_mapped_structures(), SiteCentricCorrelations(),
              AllPointCorr(), LocalCorr());

  return dict;
}

}  // namespace CASM
