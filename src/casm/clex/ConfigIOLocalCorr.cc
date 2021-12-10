#include "casm/clex/ConfigIOLocalCorr.hh"

#include "casm/app/ClexDescription.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/clex/ConfigCorrelations.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex_impl.hh"
#include "casm/clex/Supercell_impl.hh"
#include "casm/clusterography/ClusterSpecs.hh"
#include "casm/clusterography/io/json/IntegralCluster_json_io.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"

namespace CASM {
namespace ConfigIO {

// --- LocalCorr implementations -----------

const std::string LocalCorr::Name = "local_corr";

const std::string LocalCorr::Desc =
    "Local correlation values, for each equivalent phenomenal cluster, for "
    "each unit cell. The output include the unit cell indices (i,j,k), the "
    "equivalent phenomenal cluster index, phenonmenal cluster sites "
    "information, and the local correlation values. Accepts either zero "
    "arguments (`local_corr`), and uses the default cluster expansion or one "
    "argument, (`local_corr(clex_name)`) and uses the named cluster expansion. "
    "More detailed coordinates of sites for any supercell can be obtained from "
    "the `info -m NeighborListInfo` method.";

LocalCorr::LocalCorr()
    : BaseValueFormatter<jsonParser, Configuration>(Name, Desc),
      m_clex_name("") {}

jsonParser LocalCorr::evaluate(const Configuration &config) const {
  Supercell const &supercell = config.supercell();
  xtal::UnitCellCoordIndexConverter unitcellcoord_index_converter =
      supercell.sym_info().unitcellcoord_index_converter();
  std::vector<xtal::UnitCell> unitcells = xtal::make_lattice_points(
      supercell.sym_info().transformation_matrix_to_super());
  std::vector<std::vector<std::string>> site_occupant =
      xtal::allowed_molecule_names(supercell.prim());

  // result is array of:
  //
  // unitcell: [i, j, k]
  // phenomenal:
  //   sites: [[b, i, j, k], ...]
  //   linear_site_index: [l, ...]
  //   asymmetric_unit_index: [asym, ...]
  //   occupants: ["A", ...]
  //   equivalent_index: p
  // value: []

  jsonParser result = jsonParser::array();

  Index linear_unitcell_index = 0;
  for (xtal::UnitCell unitcell : unitcells) {
    Index equivalent_index = 0;
    while (equivalent_index < m_phenom_orbit.size()) {
      // get phenomenal info
      IntegralCluster clust = m_phenom_orbit[equivalent_index] + unitcell;
      std::vector<Index> asymmetric_unit_index;
      std::vector<Index> linear_site_index;
      std::vector<std::string> occupant;
      for (xtal::UnitCellCoord unitcellcoord : clust) {
        // asymmetric_unit_indices
        asymmetric_unit_index.push_back(
            m_sublattice_to_asymmetric_unit[unitcellcoord.sublattice()]);

        // linear_site_index
        Index l = unitcellcoord_index_converter(unitcellcoord);
        linear_site_index.push_back(l);

        // occupant
        Index sublattice_index = unitcellcoord.sublattice();
        Index occ = config.occ(l);
        occupant.push_back(site_occupant[sublattice_index][occ]);
      }

      jsonParser tjson;
      CASM::to_json(unitcell, tjson["unitcell"], jsonParser::as_array());
      Eigen::VectorXd corr = corr_contribution(linear_unitcell_index, config,
                                               m_clexulator[equivalent_index]);
      CASM::to_json(corr, tjson["value"], jsonParser::as_array());
      jsonParser phenom_json;
      phenom_json["sites"].put_array(clust.begin(), clust.end());
      phenom_json["asymmetric_unit_index"] = asymmetric_unit_index;
      phenom_json["linear_site_index"] = linear_site_index;
      phenom_json["occupant"] = occupant;
      phenom_json["equivalent_index"] = equivalent_index;
      tjson["phenomenal"] = phenom_json;
      result.push_back(tjson);

      ++equivalent_index;
    }
    ++linear_unitcell_index;
  }

  return result;
}

bool LocalCorr::validate(const Configuration &config) const {
  if (!m_clexulator.size()) {
    CASM::err_log() << "Error in \"local_corr\": not properly "
                       "initialized (unknown error)."
                    << std::endl;
    return false;
  }
  return true;
}

/// \brief If not yet initialized, use the default clexulator from the PrimClex
bool LocalCorr::init(const Configuration &_tmplt) const {
  if (!m_clexulator.size()) {
    const PrimClex &primclex = _tmplt.primclex();
    auto const &prim = primclex.prim();
    auto shared_prim = primclex.shared_prim();
    ClexDescription desc = m_clex_name.empty()
                               ? primclex.settings().default_clex()
                               : primclex.settings().clex(m_clex_name);
    ClexBasisSpecs const &basis_set_specs = primclex.basis_set_specs(desc.bset);
    ClusterSpecs const &cluster_specs = *basis_set_specs.cluster_specs;
    if (cluster_specs.periodicity_type() != CLUSTER_PERIODICITY_TYPE::LOCAL) {
      throw std::runtime_error(
          "Error: \"local_corr\" is only valid for local "
          "cluster expansions.");
    }

    m_clexulator = primclex.local_clexulator(desc.bset);
    m_prototype_phenom = cluster_specs.get_phenomenal_cluster();
    m_equivalents_generating_ops = make_equivalents_generating_ops(
        shared_prim, *m_prototype_phenom, cluster_specs.get_generating_group());
    m_phenom_orbit.clear();
    for (auto const &op : m_equivalents_generating_ops) {
      m_phenom_orbit.push_back(sym::copy_apply(op, *m_prototype_phenom, prim));
    }

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
/// - 'local_corr'
/// - 'local_corr(clex_name)'
bool LocalCorr::parse_args(const std::string &args) {
  std::vector<std::string> splt_vec;
  boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);

  if (!splt_vec.size()) {
    return true;
  } else if (splt_vec.size() == 1) {
    m_clex_name = splt_vec[0];
  } else {
    std::stringstream ss;
    ss << "Too many arguments for 'local_corr'.  Received: " << args << "\n";
    throw std::runtime_error(ss.str());
  }
  return true;
}

}  // namespace ConfigIO
}  // namespace CASM
