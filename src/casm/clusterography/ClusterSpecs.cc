#include "casm/clusterography/ClusterInvariants.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clusterography/IntegralCluster_impl.hh"
#include "casm/clusterography/SupercellClusterOrbits_impl.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/OrbitGeneration_impl.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/symmetry/SubOrbits_impl.hh"
#include "casm/symmetry/SupercellSymInfo.hh"

namespace CASM {

std::string ClusterSpecs::name() const { return this->_name(); }

CLUSTER_PERIODICITY_TYPE ClusterSpecs::periodicity_type() const {
  return this->_periodicity_type();
}

ClusterSpecs::PeriodicOrbitVec ClusterSpecs::make_periodic_orbits(
    IntegralClusterVec const &generating_elements) const {
  return this->_make_periodic_orbits(generating_elements);
}

ClusterSpecs::PeriodicOrbitVec ClusterSpecs::make_periodic_orbits(
    std::ostream &status) const {
  return this->_make_periodic_orbits(status);
}

ClusterSpecs::LocalOrbitVec ClusterSpecs::make_local_orbits(
    IntegralClusterVec const &generating_elements) const {
  return this->_make_local_orbits(generating_elements);
}

ClusterSpecs::LocalOrbitVec ClusterSpecs::make_local_orbits(
    std::ostream &status) const {
  return this->_make_local_orbits(status);
}

ClusterSpecs::PeriodicOrbitVec ClusterSpecs::_make_periodic_orbits(
    IntegralClusterVec const &generating_elements) const {
  throw std::runtime_error(
      "Error: make_periodic_orbits from generating elements not implemented "
      "for '" +
      name() + "'");
}

ClusterSpecs::PeriodicOrbitVec ClusterSpecs::_make_periodic_orbits(
    std::ostream &status) const {
  throw std::runtime_error("Error: make_periodic_orbits not implemented for '" +
                           name() + "'");
}

ClusterSpecs::LocalOrbitVec ClusterSpecs::_make_local_orbits(
    IntegralClusterVec const &generating_elements) const {
  throw std::runtime_error(
      "Error: make_local_orbits from generating elements not implemented for "
      "'" +
      name() + "'");
}

ClusterSpecs::LocalOrbitVec ClusterSpecs::_make_local_orbits(
    std::ostream &status) const {
  throw std::runtime_error("Error: make_local_orbits not implemented for '" +
                           name() + "'");
}

void ClusterSpecs::set_generating_group(SymGroup const &_generating_group) {
  return this->_set_generating_group(_generating_group);
}

SymGroup const &ClusterSpecs::get_generating_group() const {
  return this->_get_generating_group();
}

IntegralCluster const &ClusterSpecs::get_phenomenal_cluster() const {
  return this->_get_phenomenal_cluster();
}

void ClusterSpecs::_set_generating_group(SymGroup const &_generating_group) {
  throw std::runtime_error("Error: set_generating_group not implemented for '" +
                           name() + "'");
}

SymGroup const &ClusterSpecs::_get_generating_group() const {
  throw std::runtime_error("Error: get_generating_group not implemented for '" +
                           name() + "'");
}

IntegralCluster const &ClusterSpecs::_get_phenomenal_cluster() const {
  throw std::runtime_error(
      "Error: get_phenomenal_cluster not implemented for '" + name() + "'");
}

const std::string PeriodicMaxLengthClusterSpecs::method_name =
    "periodic_max_length";

PeriodicMaxLengthClusterSpecs::PeriodicMaxLengthClusterSpecs(
    std::shared_ptr<Structure const> _shared_prim,
    SymGroup const &_generating_group, SiteFilterFunction const &_site_filter,
    std::vector<double> const &_max_length,
    std::vector<IntegralClusterOrbitGenerator> const &_custom_generators)
    : shared_prim(_shared_prim),
      generating_group(_generating_group),
      sym_compare(shared_prim, shared_prim->lattice().tol()),
      site_filter(_site_filter),
      max_length(_max_length),
      custom_generators(_custom_generators) {
  if (max_length.size() == 0) {
    throw libcasm_runtime_error(
        "Error in PeriodicMaxLengthClusterSpecs: max_length.size() == 0 (must "
        "be greater than 0).");
  }
}

std::string PeriodicMaxLengthClusterSpecs::_name() const { return method_name; }

CLUSTER_PERIODICITY_TYPE PeriodicMaxLengthClusterSpecs::_periodicity_type()
    const {
  return CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC;
};

ClusterSpecs::PeriodicOrbitVec
PeriodicMaxLengthClusterSpecs::_make_periodic_orbits(
    IntegralClusterVec const &generating_elements) const {
  return generate_orbits(generating_elements, generating_group, sym_compare);
}

ClusterSpecs::PeriodicOrbitVec
PeriodicMaxLengthClusterSpecs::_make_periodic_orbits(
    std::ostream &status) const {
  typedef PrimPeriodicOrbit<IntegralCluster> orbit_type;
  std::vector<OrbitBranchSpecs<orbit_type> > specs;

  for (int branch = 0; branch < max_length.size(); ++branch) {
    std::vector<xtal::UnitCellCoord> candidate_sites;
    CandidateSitesFunction f;
    if (branch == 0) {
      f = empty_neighborhood();
    } else if (branch == 1) {
      f = origin_neighborhood();
    } else {
      f = max_length_neighborhood(max_length[branch]);
    }
    candidate_sites = f(*shared_prim, site_filter);

    ClusterFilterFunction cluster_filter;
    if (branch <= 1) {
      cluster_filter = all_clusters_filter();
    } else {
      cluster_filter = max_length_cluster_filter(max_length[branch]);
    }

    specs.emplace_back(*shared_prim, candidate_sites.begin(),
                       candidate_sites.end(), generating_group, cluster_filter,
                       sym_compare);
  }

  // now generate orbits
  PeriodicOrbitVec orbits;
  make_orbits(specs.begin(), specs.end(), custom_generators,
              std::back_inserter(orbits), status);
  return orbits;
}

void PeriodicMaxLengthClusterSpecs::_set_generating_group(
    SymGroup const &_generating_group) {
  this->generating_group = _generating_group;
}

SymGroup const &PeriodicMaxLengthClusterSpecs::_get_generating_group() const {
  return generating_group;
}

std::string const LocalMaxLengthClusterSpecs::method_name = "local_max_length";

LocalMaxLengthClusterSpecs::LocalMaxLengthClusterSpecs(
    std::shared_ptr<Structure const> _shared_prim,
    SymGroup const &_generating_group, IntegralCluster const &_phenomenal,
    SiteFilterFunction const &_site_filter,
    std::vector<double> const &_max_length,
    std::vector<double> const &_cutoff_radius, bool _include_phenomenal_sites,
    std::vector<IntegralClusterOrbitGenerator> const &_custom_generators)
    : shared_prim(_shared_prim),
      generating_group(_generating_group),
      sym_compare(shared_prim, shared_prim->lattice().tol()),
      phenomenal(_phenomenal),
      site_filter(_site_filter),
      max_length(_max_length),
      cutoff_radius(_cutoff_radius),
      include_phenomenal_sites(_include_phenomenal_sites),
      custom_generators(_custom_generators) {
  if (max_length.size() == 0) {
    throw libcasm_runtime_error(
        "Error in LocalMaxLengthClusterSpecs: max_length.size() == 0 (must be "
        "greater than 0).");
  }
  if (max_length.size() != cutoff_radius.size()) {
    throw libcasm_runtime_error(
        "Error in LocalMaxLengthClusterSpecs: max_length.size() != "
        "cutoff_radius.size() (must be equal).");
  }
}

std::string LocalMaxLengthClusterSpecs::_name() const { return method_name; }

CLUSTER_PERIODICITY_TYPE LocalMaxLengthClusterSpecs::_periodicity_type() const {
  return CLUSTER_PERIODICITY_TYPE::LOCAL;
};

ClusterSpecs::LocalOrbitVec LocalMaxLengthClusterSpecs::_make_local_orbits(
    IntegralClusterVec const &generating_elements) const {
  return generate_orbits(generating_elements, generating_group, sym_compare);
}

ClusterSpecs::LocalOrbitVec LocalMaxLengthClusterSpecs::_make_local_orbits(
    std::ostream &status) const {
  typedef LocalOrbit<IntegralCluster> orbit_type;
  std::vector<OrbitBranchSpecs<orbit_type> > specs;

  for (int branch = 0; branch < max_length.size(); ++branch) {
    std::vector<xtal::UnitCellCoord> candidate_sites;
    CandidateSitesFunction f;
    if (branch == 0) {
      f = empty_neighborhood();
    } else {
      f = cutoff_radius_neighborhood(phenomenal, cutoff_radius[branch],
                                     include_phenomenal_sites);
    }
    candidate_sites = f(*shared_prim, site_filter);

    ClusterFilterFunction cluster_filter;
    if (branch <= 1) {
      cluster_filter = all_clusters_filter();
    } else {
      cluster_filter = max_length_cluster_filter(max_length[branch]);
    }

    specs.emplace_back(*shared_prim, candidate_sites.begin(),
                       candidate_sites.end(), generating_group, cluster_filter,
                       sym_compare);
  }

  // now generate orbits
  LocalOrbitVec orbits;
  make_orbits(specs.begin(), specs.end(), custom_generators,
              std::back_inserter(orbits), status);
  return orbits;
}

void LocalMaxLengthClusterSpecs::_set_generating_group(
    SymGroup const &_generating_group) {
  this->generating_group = _generating_group;
}

SymGroup const &LocalMaxLengthClusterSpecs::_get_generating_group() const {
  return generating_group;
}

IntegralCluster const &LocalMaxLengthClusterSpecs::_get_phenomenal_cluster()
    const {
  return phenomenal;
}

/// \brief A set of SymOp that construct distinct equivalent phenomenal clusters
///
/// \param shared_prim Prim structure
/// \param phenomenal Phenomenal cluster about which local clusters are
///     generated
/// \param generating_group Local clusters generating group. Must leave the
///     phenomenal cluster sites invariant.
std::vector<SymOp> make_equivalents_generating_ops(
    std::shared_ptr<Structure const> const &shared_prim,
    IntegralCluster const &phenomenal, SymGroup const &generating_group) {
  PrimPeriodicSymCompare<IntegralCluster> periodic_sym_compare(
      shared_prim, shared_prim->lattice().tol());

  // Invariant group of phenomenal cluster sites, excluding any other
  // considerations like symmetry-breaking due to hop type
  SymGroup sites_invariant_group = make_invariant_subgroup(
      phenomenal, shared_prim->factor_group(), periodic_sym_compare);

  // If generating_group is a subgroup of sites_invariant_group, then there
  // may be some distinct phenomenal clusters on the same sites (ex: distinct
  // orientations of local cluster expansions about the same phenomenal cluster)
  std::set<SymOp, SubOrbits_impl::OpCompare> equivalents_on_phenomenal;
  for (auto const &op : sites_invariant_group) {
    Index min = op.index();
    SymOp min_op = op;
    for (auto const &subgroup_op : generating_group) {
      SymOp product = subgroup_op * op;
      if (product.index() < min_op.index()) {
        min = product.index();
        min_op = product;
      }
    }
    equivalents_on_phenomenal.insert(min_op);
  }

  // generate the phenomenal cluster orbit
  PrimPeriodicOrbit<IntegralCluster> phenomenal_cluster_orbit(
      phenomenal, shared_prim->factor_group(), periodic_sym_compare);
  auto const &equivalence_map = phenomenal_cluster_orbit.equivalence_map();

  // find the phenomenal cluster in the orbit,
  // get the operation that transforms phenomenal to prototype
  Index phenomenal_index = -1;
  IntegralCluster sorted_phenom =
      periodic_sym_compare.representation_prepare(phenomenal);
  for (int i = 0; i < equivalence_map.size(); ++i) {
    IntegralCluster const &equiv = phenomenal_cluster_orbit.element(i);
    IntegralCluster sorted_equiv =
        periodic_sym_compare.representation_prepare(equiv);
    if (periodic_sym_compare.equal(sorted_phenom, sorted_equiv)) {
      phenomenal_index = i;
    }
  }
  if (phenomenal_index == -1) {
    throw std::runtime_error(
        "Error constructing a local ClexBasis: the phenomenal "
        "cluster must be in the origin unit cell orbit.");
  }

  SymOp to_prototype_op = equivalence_map[phenomenal_index][0].inverse();

  // operations to generate all equivalent clexulators at all clusters in
  // phenomenal cluster orbit
  std::set<SymOp, SubOrbits_impl::OpCompare> equivalents_generating_ops;
  for (SymOp const &cluster_op : equivalents_on_phenomenal) {
    for (int equivalent_index = 0; equivalent_index < equivalence_map.size();
         ++equivalent_index) {
      equivalents_generating_ops.insert(equivalence_map[equivalent_index][0] *
                                        to_prototype_op * cluster_op);
    }
  }

  return std::vector<SymOp>(equivalents_generating_ops.begin(),
                            equivalents_generating_ops.end());
}

/// \brief The extended equivalence map maps a prototype to equivalent
///     clusters in the current and equivalent orbitrees
///
/// Usage:
/// \code
/// IntegralCluster equiv_cluster = sym_compare.copy_apply(
///     extended_equivalence_map[phenom_index][cluster_index][symop_index],
///     prototype)
/// \endcode
/// Where `equiv_cluster` is the cluster_index-th cluster around the
/// phenom_index-th equivalent phenomenal cluster, for all symop_index.
///
multivector<SymOp>::X<3> make_extended_equivalence_map(
    multivector<SymOp>::X<2> const &equivalence_map,
    std::vector<SymOp> const &equivalents_generating_ops) {
  multivector<SymOp>::X<3> extended_equivalence_map;
  for (SymOp const &equiv_phenom_op : equivalents_generating_ops) {
    multivector<SymOp>::X<2> equiv_equivalence_map;
    // now use equiv_local_op and the local cluster orbit equivalence map to
    // generate new rows in the exteneded equivalence map corresponding to local
    // clusters around `equiv_phenomenal`
    for (auto const &row : equivalence_map) {
      std::vector<SymOp> extended_equivalence_map_row;
      for (SymOp const &local_orbit_op : row) {
        extended_equivalence_map_row.push_back(equiv_phenom_op *
                                               local_orbit_op);
      }
      equiv_equivalence_map.push_back(extended_equivalence_map_row);
    }
    extended_equivalence_map.push_back(equiv_equivalence_map);
  }
  return extended_equivalence_map;
}

GenericPeriodicClusterSpecs::GenericPeriodicClusterSpecs(
    std::string _method_name, std::shared_ptr<Structure const> _shared_prim,
    SymGroup const &_generating_group, SymCompareType const &_sym_compare,
    SiteFilterFunction _site_filter,
    std::vector<ClusterFilterFunction> _cluster_filter,
    std::vector<CandidateSitesFunction> _candidate_sites,
    std::vector<IntegralClusterOrbitGenerator> _custom_generators)
    : shared_prim(_shared_prim),
      generating_group(_generating_group),
      sym_compare(_sym_compare),
      site_filter(_site_filter),
      cluster_filter(_cluster_filter),
      candidate_sites(_candidate_sites),
      custom_generators(_custom_generators),
      m_method_name(_method_name) {
  if (candidate_sites.size() == 0) {
    throw libcasm_runtime_error(
        "Error in GenericPeriodicClusterSpecs: candidate_sites.size() == 0 "
        "(must be greater than 0).");
  }
}

std::string GenericPeriodicClusterSpecs::_name() const { return m_method_name; }

CLUSTER_PERIODICITY_TYPE GenericPeriodicClusterSpecs::_periodicity_type()
    const {
  return CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC;
}

ClusterSpecs::PeriodicOrbitVec
GenericPeriodicClusterSpecs::_make_periodic_orbits(
    IntegralClusterVec const &generating_elements) const {
  return generate_orbits(generating_elements, generating_group, sym_compare);
}

ClusterSpecs::PeriodicOrbitVec
GenericPeriodicClusterSpecs::_make_periodic_orbits(std::ostream &status) const {
  if (cluster_filter.size() != candidate_sites.size()) {
    throw std::runtime_error(
        "Error in GenericPeriodicClusterSpecs::_make_periodic_orbits: "
        "cluster_filter.size() != candidate_sites.size()");
  }
  typedef PrimPeriodicOrbit<IntegralCluster> OrbitType;
  std::vector<OrbitBranchSpecs<OrbitType> > specs;

  for (int branch = 0; branch < cluster_filter.size(); ++branch) {
    std::vector<xtal::UnitCellCoord> tmp;
    tmp = candidate_sites[branch](*shared_prim, site_filter);

    specs.emplace_back(*shared_prim, tmp.begin(), tmp.end(), generating_group,
                       cluster_filter[branch], sym_compare);
  }

  // now generate orbits
  PeriodicOrbitVec orbits;
  make_orbits(specs.begin(), specs.end(), custom_generators,
              std::back_inserter(orbits), status);
  return orbits;
}

void GenericPeriodicClusterSpecs::_set_generating_group(
    SymGroup const &_generating_group) {
  this->generating_group = _generating_group;
}

SymGroup const &GenericPeriodicClusterSpecs::_get_generating_group() const {
  return generating_group;
}

GenericLocalClusterSpecs::GenericLocalClusterSpecs(
    std::string _method_name, std::shared_ptr<Structure const> _shared_prim,
    SymGroup const &_generating_group, SymCompareType const &_sym_compare,
    SiteFilterFunction _site_filter,
    std::vector<ClusterFilterFunction> _cluster_filter,
    std::vector<CandidateSitesFunction> _candidate_sites,
    std::vector<IntegralClusterOrbitGenerator> _custom_generators)
    : shared_prim(_shared_prim),
      generating_group(_generating_group),
      sym_compare(_sym_compare),
      site_filter(_site_filter),
      cluster_filter(_cluster_filter),
      candidate_sites(_candidate_sites),
      custom_generators(_custom_generators),
      m_method_name(_method_name) {
  if (candidate_sites.size() == 0) {
    throw libcasm_runtime_error(
        "Error in GenericLocalClusterSpecs: candidate_sites.size() == 0 (must "
        "be greater than 0).");
  }
  if (candidate_sites.size() != cluster_filter.size()) {
    throw libcasm_runtime_error(
        "Error in GenericLocalClusterSpecs: candidate_sites.size() != "
        "cluster_filter.size() (must be equal).");
  }
}

std::string GenericLocalClusterSpecs::_name() const { return m_method_name; }

CLUSTER_PERIODICITY_TYPE GenericLocalClusterSpecs::_periodicity_type() const {
  return CLUSTER_PERIODICITY_TYPE::LOCAL;
}

ClusterSpecs::LocalOrbitVec GenericLocalClusterSpecs::_make_local_orbits(
    IntegralClusterVec const &generating_elements) const {
  return generate_orbits(generating_elements, generating_group, sym_compare);
}

ClusterSpecs::LocalOrbitVec GenericLocalClusterSpecs::_make_local_orbits(
    std::ostream &status) const {
  if (cluster_filter.size() != candidate_sites.size()) {
    throw std::runtime_error(
        "Error in GenericLocalClusterSpecs::_make_local_orbits: "
        "candidate_sites.size() != cluster_filter.size()");
  }
  typedef LocalOrbit<IntegralCluster> OrbitType;
  std::vector<OrbitBranchSpecs<OrbitType> > specs;

  for (int branch = 0; branch < cluster_filter.size(); ++branch) {
    std::vector<xtal::UnitCellCoord> tmp;
    tmp = candidate_sites[branch](*shared_prim, site_filter);

    specs.emplace_back(*shared_prim, tmp.begin(), tmp.end(), generating_group,
                       cluster_filter[branch], sym_compare);
  }

  // now generate orbits
  LocalOrbitVec orbits;
  make_orbits(specs.begin(), specs.end(), custom_generators,
              std::back_inserter(orbits), status);
  return orbits;
}

void GenericLocalClusterSpecs::_set_generating_group(
    SymGroup const &_generating_group) {
  this->generating_group = _generating_group;
}

SymGroup const &GenericLocalClusterSpecs::_get_generating_group() const {
  return generating_group;
}

namespace ClusterSpecs_impl {

class DoFSitesFilter {
 public:
  DoFSitesFilter(std::vector<DoFKey> const &_dofs) : dofs(_dofs) {}

  bool operator()(xtal::Site const &site) {
    if (dofs.empty() &&
        (site.dof_size() != 0 || site.occupant_dof().size() > 1)) {
      return true;
    }
    for (DoFKey const &dof : dofs) {
      if (site.has_dof(dof)) {
        return true;
      } else if (dof == "occ" && site.occupant_dof().size() > 1) {
        return true;
      }
    }
    return false;
  }

  std::vector<DoFKey> dofs;
};

class AllClusters {
 public:
  bool operator()(IntegralCluster const &clust) { return true; }
};

class MaxLengthClusterFilter {
 public:
  MaxLengthClusterFilter(double _max_length) : max_length(_max_length) {}

  bool operator()(IntegralCluster const &clust) {
    if (clust.size() <= 1) {
      return true;
    }
    ClusterInvariants invariants{clust};
    return invariants.displacement().back() < max_length;
  }

 private:
  double max_length;
};

class EmptyNeighborhood {
 public:
  std::vector<xtal::UnitCellCoord> operator()(Structure const &prim,
                                              SiteFilterFunction site_filter) {
    return std::vector<xtal::UnitCellCoord>{};
  }
};

class OriginNeighborhood {
 public:
  std::vector<xtal::UnitCellCoord> operator()(Structure const &prim,
                                              SiteFilterFunction site_filter) {
    std::vector<xtal::UnitCellCoord> result;
    for (int i = 0; i < prim.basis().size(); ++i) {
      if (site_filter(prim.basis()[i])) {
        result.emplace_back(i, 0, 0, 0);
      }
    }
    return result;
  }
};

class MaxLengthNeighborhood {
 public:
  MaxLengthNeighborhood(double _max_length) : max_length(_max_length){};

  std::vector<xtal::UnitCellCoord> operator()(Structure const &prim,
                                              SiteFilterFunction site_filter) {
    std::vector<xtal::UnitCellCoord> result;
    double xtal_tol = prim.lattice().tol();
    neighborhood(prim, max_length, site_filter, std::back_inserter(result),
                 xtal_tol);
    return result;
  }

 private:
  double max_length;
};

// class ScelNeighborhood {
// public:
//   ScelNeighborhood(Eigen::Matrix3l const &supercell_matrix):
//     lattice_points(xtal::make_lattice_points(supercell_matrix)) {}
//
//   std::vector<xtal::UnitCellCoord> operator()(Structure const &prim,
//   SiteFilterFunction site_filter) {
//     std::vector<xtal::UnitCellCoord> result;
//     int b = 0;
//     for(auto const &site : prim.basis()) {
//       if(site_filter(site)) {
//         for(auto const &lattice_point : lattice_points) {
//           result.emplace_back(b, lattice_point);
//         }
//       }
//       ++b;
//     }
//     return result;
//   }
//
// private:
//   std::vector<xtal::UnitCell> lattice_points;
// };

/// Generate a vector of UnitCellCoord that are within cutoff_radius distance to
/// any site in the phenomenal cluster
class CutoffRadiusNeighborhood {
 public:
  CutoffRadiusNeighborhood(IntegralCluster const &_phenomenal,
                           double _cutoff_radius,
                           bool _include_phenomenal_sites)
      : phenomenal(_phenomenal),
        cutoff_radius(_cutoff_radius),
        include_phenomenal_sites(_include_phenomenal_sites){};

  std::vector<xtal::UnitCellCoord> operator()(Structure const &prim,
                                              SiteFilterFunction site_filter) {
    std::vector<xtal::UnitCellCoord> result;
    double xtal_tol = prim.lattice().tol();
    neighborhood(phenomenal, cutoff_radius, site_filter,
                 include_phenomenal_sites, std::back_inserter(result),
                 xtal_tol);
    return result;
  }

 private:
  IntegralCluster phenomenal;
  double cutoff_radius;
  bool include_phenomenal_sites;
};

}  // namespace ClusterSpecs_impl

/// \brief Generate clusters using all Site
bool all_sites_filter(xtal::Site const &site) { return true; }

/// \brief Generate clusters using Site with site_occupant.size() > 1
bool alloy_sites_filter(xtal::Site const &site) {
  return site.occupant_dof().size() > 1;
}

/// \brief Generate clusters using Site with specified DoF
///
/// If dofs is empty, return true if Site has any continuous DoF or >1 allowed
/// occupant DoF If dofs is not empty, return true if Site has any of the DoF
/// types included. Use "occ" for / Site with >1 occupant allowed
SiteFilterFunction dof_sites_filter(std::vector<DoFKey> const &dofs) {
  return ClusterSpecs_impl::DoFSitesFilter{dofs};
}

/// Accept all clusters
ClusterFilterFunction all_clusters_filter() {
  return ClusterSpecs_impl::AllClusters{};
}

/// Accept clusters with max pair distance less than max_length
ClusterFilterFunction max_length_cluster_filter(double max_length) {
  return ClusterSpecs_impl::MaxLengthClusterFilter{max_length};
}

/// No sites (for null orbit, or global dof only)
CandidateSitesFunction empty_neighborhood() {
  return ClusterSpecs_impl::EmptyNeighborhood{};
}

/// Only sites in the origin unit cell {b, 0, 0, 0}
CandidateSitesFunction origin_neighborhood() {
  return ClusterSpecs_impl::OriginNeighborhood{};
}

/// Sites within max_length distance to any site in the origin unit cell {b, 0,
/// 0, 0}
CandidateSitesFunction max_length_neighborhood(double max_length) {
  return ClusterSpecs_impl::MaxLengthNeighborhood{max_length};
}

/// Sites within cutoff_radius distance to any site in the phenomenal cluster
CandidateSitesFunction cutoff_radius_neighborhood(
    IntegralCluster const &phenomenal, double cutoff_radius,
    bool include_phenomenal_sites) {
  return ClusterSpecs_impl::CutoffRadiusNeighborhood{phenomenal, cutoff_radius,
                                                     include_phenomenal_sites};
}

}  // namespace CASM
