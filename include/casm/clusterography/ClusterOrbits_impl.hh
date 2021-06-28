#ifndef CASM_ClusterOrbits_impl
#define CASM_ClusterOrbits_impl

#include "casm/clusterography/ClusterInvariants.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/clusterography/ClusterSymCompare.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clusterography/IntegralClusterSymCompareTraits_impl.hh"
#include "casm/clusterography/SubClusterGenerator.hh"
#include "casm/container/Counter.hh"
#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/OrbitGeneration_impl.hh"
#include "casm/symmetry/Orbit_impl.hh"

namespace CASM {

/* -- Cluster Orbit generating function definitions
 * ------------------------------------- */

/// \brief Output the neighborhood of UnitCellCoord within max_radius of any
/// site in unit cell
///
/// \param unit The unit cell Structure
/// \param max_radius The neighborhood distance cutoff
/// \param site_filter A filter function that returns true for CoordType that
///        should be considered for the neighborhood
/// \param result Output iterator for container of UnitCellCoord
/// \param xtal_tol Crystallography tolerance used to contstruct UnitCellCoord
/// from CoordType
///
/// \returns Output iterator after generating the neighborhood
///
/// \ingroup IntegralCluster
///
template <typename OutputIterator>
OutputIterator neighborhood(Structure const &unit, double max_radius,
                            SiteFilterFunction site_filter,
                            OutputIterator result, double xtal_tol) {
  auto dim = unit.lattice().enclose_sphere(max_radius);
  EigenCounter<Eigen::Vector3i> grid_count(-dim, dim,
                                           Eigen::Vector3i::Constant(1));
  xtal::Coordinate lat_point(unit.lattice());
  const auto &basis = unit.basis();

  do {
    lat_point.frac() = grid_count().cast<double>();

    for (auto it = basis.begin(); it != basis.end(); ++it) {
      if (!site_filter(*it)) {
        continue;
      }

      xtal::Coordinate test(*it + lat_point);
      auto within_radius = [&](const xtal::Coordinate &coord) {
        return test.dist(coord) < max_radius;
      };
      if (std::any_of(basis.begin(), basis.end(), within_radius)) {
        *result++ = xtal::UnitCellCoord::from_coordinate(unit, test, xtal_tol);
      }
    }
  } while (++grid_count);
  return result;
}

/// \brief Output the neighborhood of sites within cutoff_radius of any sites in
/// the phenomenal
///
/// \param phenomenal IntegralCluster
/// \param cutoff_radius The neighborhood distance cutoff
/// \param site_filter A filter function that returns true for UnitCellCoord
/// that
///        should be considered for the neighborhood
/// \param result Output iterator for container of UnitCellCoord
/// \param xtal_tol Crystallography tolerance used to contstruct UnitCellCoord
///
/// \returns Output iterator after generating the neighborhood
///
/// \ingroup IntegralCluster
///
template <typename OutputIterator>
OutputIterator neighborhood(IntegralCluster const &phenomenal,
                            double cutoff_radius,
                            SiteFilterFunction site_filter,
                            bool include_phenomenal_sites,
                            OutputIterator result, double xtal_tol) {
  const auto &prim = phenomenal.prim();

  int max_low_shift = 0;
  int max_high_shift = 0;
  for (auto it = phenomenal.begin(); it != phenomenal.end(); it++) {
    Eigen::Vector3l vec = it->unitcell();
    if (vec.maxCoeff() > max_high_shift) {
      max_high_shift = vec.maxCoeff();
    }
    if (vec.minCoeff() < max_low_shift) {
      max_low_shift = vec.minCoeff();
    }
  }

  // make grid counter
  auto dim = prim.lattice().enclose_sphere(cutoff_radius);
  Eigen::Vector3i ones(1, 1, 1);
  EigenCounter<Eigen::Vector3i> grid_count(
      -dim + (max_low_shift * ones).cast<int>(),
      dim + (max_high_shift * ones).cast<int>(), Eigen::Vector3i::Constant(1));

  /// lattice scaling
  xtal::Coordinate lat_point(prim.lattice());

  const auto &basis = prim.basis();

  do {
    lat_point.frac() = grid_count().cast<double>();

    for (auto it = basis.begin(); it != basis.end(); ++it) {
      if (!site_filter(*it)) {
        continue;
      }

      xtal::Coordinate test{*it + lat_point};

      if (!include_phenomenal_sites) {
        auto is_almost_equal = [&](const xtal::UnitCellCoord &uccoord) {
          return test.dist(uccoord.coordinate(prim)) < xtal_tol;
        };
        if (std::any_of(phenomenal.begin(), phenomenal.end(),
                        is_almost_equal)) {
          continue;
        }
      }

      auto within_radius = [&](const xtal::UnitCellCoord &uccoord) {
        return test.dist(uccoord.coordinate(prim)) < cutoff_radius;
      };
      if (std::any_of(phenomenal.begin(), phenomenal.end(), within_radius)) {
        *result++ = xtal::UnitCellCoord::from_coordinate(prim, test, xtal_tol);
      }
    }
  } while (++grid_count);
  return result;
}

/// \brief Iterate over all sites in an orbit and insert a UnitCellCoord
///
/// \param orbit an Orbit<IntegralCluster>
/// \param result an OutputIterator for UnitCellCoord
///
/// \result the resulting OutputIterator
///
/// This simply outputs all UnitCellCoord in all equivalent clusters
///
/// \ingroup IntegralCluster
///
template <typename OrbitType, typename OutputIterator>
OutputIterator local_orbit_neighborhood(const OrbitType &orbit,
                                        OutputIterator result) {
  for (const auto &equiv : orbit) {
    for (const auto &site : equiv) {
      *result++ = site;
    }
  }
  return result;
}

/// \brief Iterate over all sites in all orbits and insert a UnitCellCoord
///
/// \param begin,end Range of Orbit<IntegralCluster>
/// \param result an OutputIterator for UnitCellCoord
///
/// This simply outputs all UnitCellCoord in all equivalent clusters of each
/// orbit
///
/// \ingroup IntegralCluster
///
template <typename ClusterOrbitIterator, typename OutputIterator>
OutputIterator local_neighborhood(ClusterOrbitIterator begin,
                                  ClusterOrbitIterator end,
                                  OutputIterator result) {
  // create a neighborhood of all UnitCellCoord that an Orbitree touches
  for (auto it = begin; it != end; ++it) {
    result = local_orbit_neighborhood(*it, result);
  }
  return result;
}

/// \brief Check if periodic images of sites in an orbit are overlapping in
/// supercells defined by the given superlattice transformation matrix
///
/// \param local_orbits Vector of IntegralCluster orbits defining the local
/// neighborhood \param transf_mat A supercell transformation matrix
///
template <typename OrbitType>
bool has_local_neighborhood_overlap(std::vector<OrbitType> const &local_orbits,
                                    Eigen::Matrix3i const &transf_mat) {
  // use set of pair: <b, lattice site index> to check for site overlap
  std::set<std::pair<Index, Index>> present;
  xtal::LinearIndexConverter<xtal::UnitCell> convert{transf_mat};
  for (auto const &orbit : local_orbits) {
    for (auto const &cluster : orbit) {
      for (auto const &uccoord : cluster) {
        auto res =
            present.emplace(uccoord.sublattice(), convert[uccoord.unitcell()]);
        // if no insert took place, there is an overlap, so return true
        if (!res.second) {
          return true;
        }
      }
    }
  }
  // if no set insertion collision then no overlap
  return false;
}

/// Make site dependency neighborhoods
///
/// \param begin Iterator at the beginning of a cluster orbit range
/// \param end Iterator at the end of a cluster orbit range
///
/// \returns a map of UnitCellCoord (translated as necessary to the canonical
/// unit cell by the SymCompare method), to the set of UnitCellCoord that it is
/// in clusters with. For cluster function evaluation, this gives the
/// neighborhood of sites whose DoF values are needed to evaluate clusters
/// involving a particular site.
///
/// Note:
/// - keys of result are guaranteed to be in canonical translation unit
template <typename ClusterOrbitIterator>
std::map<xtal::UnitCellCoord, std::set<xtal::UnitCellCoord>>
make_site_dependency_neighborhoods(ClusterOrbitIterator begin,
                                   ClusterOrbitIterator end) {
  std::map<xtal::UnitCellCoord, std::set<xtal::UnitCellCoord>> result;

  if (begin == end) return result;

  typedef IntegralCluster cluster_type;
  typedef typename ClusterOrbitIterator::value_type orbit_type;

  xtal::UnitCellCoord const *ucc_ptr(nullptr);
  ClusterOrbitIterator begin2 = begin;
  for (; begin2 != end; ++begin2) {
    auto const &orbit = *begin2;
    for (auto const &equiv : orbit) {
      for (xtal::UnitCellCoord const &ucc : equiv) {
        ucc_ptr = &ucc;
        break;
      }
      if (ucc_ptr) break;
    }
    if (ucc_ptr) break;
  }
  if (!ucc_ptr) return result;

  Structure const &prim(begin->prototype().prim());
  SymGroup identity_group(prim.factor_group().begin(),
                          (prim.factor_group().begin()) + 1);
  orbit_type empty_orbit(cluster_type(prim), identity_group,
                         begin->sym_compare());

  // Loop over each site in each cluster of each orbit
  for (; begin != end; ++begin) {
    auto const &orbit = *begin;
    for (auto const &equiv : orbit) {
      for (xtal::UnitCellCoord const &ucc : equiv) {
        // create a test cluster from prototype
        cluster_type test(empty_orbit.prototype());

        // add the new site
        test.elements().push_back(ucc);
        test = orbit.sym_compare().prepare(test);

        xtal::UnitCell trans = test.element(0).unitcell() - ucc.unitcell();
        for (xtal::UnitCellCoord const &ucc2 : equiv) {
          result[test.element(0)].insert(ucc2 + trans);
        }
      }
    }
  }
  return result;
}

/// \brief Return superlattice transf. matrices for which
/// has_local_neighborhood_overlap is false
///
/// \param local_orbits Vector of IntegralCluster orbits defining the local
/// neighborhood \param transf_mat_options A vector of supercell transformation
/// matrices
///
template <typename OrbitType>
std::vector<Eigen::Matrix3i> viable_supercells(
    std::vector<OrbitType> &local_orbits,
    const std::vector<Eigen::Matrix3i> &transf_mat_options) {
  std::vector<Eigen::Matrix3i> results;
  for (auto &transf_mat : transf_mat_options) {
    if (!has_local_neighborhood_overlap(local_orbits, transf_mat)) {
      results.push_back(transf_mat);
    }
  }
  return results;
}

/* -- Custom clusters --- */

/// \brief Given a cluster, generate all subcluster generators
///
/// \param cluster A cluster to generate subclusters of
/// \param generators An OrbitGeneratorSet<OrbitType>& to store generating
///        elements for subclusters of cluster
/// \param stutus Stream for status messages
///
/// Uses SymCompareType::compare to find unique generating elements
///
/// \ingroup IntegralCluster
///
template <typename OrbitType>
OrbitGenerators<OrbitType> &insert_subcluster_generators(
    typename OrbitType::Element cluster,
    OrbitGenerators<OrbitType> &generators) {
  typedef typename OrbitType::Element cluster_type;
  typedef OrbitType orbit_type;

  // Construct a functor that returns takes a cluster and returns it in a
  // canonical form
  CanonicalGenerator<orbit_type> canonical_generator(generators.group,
                                                     generators.sym_compare);

  // SubClusterGenerator allows iterating over subclusters (includes null and
  // original cluster)
  SubClusterGenerator<cluster_type> it(cluster);
  SubClusterGenerator<cluster_type> end;

  while (it != end) {
    generators.elements.insert(canonical_generator(*it++));
  }

  return generators;
}

/* -- Cluster Orbit access/usage function definitions
 * ------------------------------------- */

/// \brief Returns the first range containing orbits of the requested orbit
/// branch in the given range of Orbit
///
/// \ingroup Clusterography
///
template <typename OrbitIterator>
std::pair<OrbitIterator, OrbitIterator> orbit_branch(OrbitIterator begin,
                                                     OrbitIterator end,
                                                     Index size) {
  auto branch_begin = std::find_if(
      begin, end, [&](const typename OrbitIterator::value_type &orbit) {
        return orbit.prototype().size() == size;
      });

  auto branch_end = std::find_if(
      branch_begin, end, [&](const typename OrbitIterator::value_type &orbit) {
        return orbit.prototype().size() == size + 1;
      });

  return std::pair<OrbitIterator, OrbitIterator>(branch_begin, branch_end);
}

// ------- Generating elements generation ------------------------------------

namespace {

/// \brief Insert the null cluster
///
/// \param prim A PrimType
/// \param generators OrbitGenerators<OrbitType> instance collecting orbit
///        generating elements
///
/// \relates IntegralCluster
///
template <typename OrbitType>
OrbitGenerators<OrbitType> &_insert_null_cluster_generator(
    const Structure &prim, OrbitGenerators<OrbitType> &generators) {
  typedef typename OrbitType::Element cluster_type;

  // create a prototype cluster
  cluster_type test(prim);
  generators.insert(test);

  return generators;
}

/// \brief Insert the generating elements for the asymmetric unit, including all
/// sites
///
/// \param prim A PrimType
/// \param generators OrbitGenerators<OrbitType> instance collecting orbit
///        generating elements
///
/// \relates IntegralCluster
///
template <typename OrbitType>
OrbitGenerators<OrbitType> &_insert_asymmetric_unit_generators(
    const Structure &prim, OrbitGenerators<OrbitType> &generators) {
  typedef typename OrbitType::Element cluster_type;

  // for all sites in the basis
  for (int i = 0; i < prim.basis().size(); i++) {
    // create a prototype cluster
    cluster_type test(prim);
    test.elements().emplace_back(i, xtal::UnitCell(0, 0, 0));
    generators.insert(test);
  }

  return generators;
}

/// \brief Use orbits of size n to insert generating elements for orbits of size
/// n+1
///
/// \param begin,end A range of input orbit generating clusters of size n
/// \param specs OrbitBranchSpecs for orbits of size n+1
/// \param generators An OrbitGeneratorSet<OrbitType> >& to store
///        generating elements for orbits of size n+1
/// \param stutus Stream for status messages
///
/// Uses SymCompareType::compare to find unique generating elements
///
/// \ingroup IntegralCluster
///
template <typename OrbitType, typename OrbitGeneratorIterator>
OrbitGenerators<OrbitType> &_insert_next_orbitbranch_generators(
    OrbitGeneratorIterator begin, OrbitGeneratorIterator end,
    const OrbitBranchSpecs<OrbitType> &specs,
    OrbitGenerators<OrbitType> &generators, std::ostream &status) {
  typedef typename OrbitType::Element cluster_type;

  const auto &filter = specs.filter();

  // print status messages
  std::string clean(100, ' ');

  // contains a pair of iterators over candidate UnitCellCoord
  auto candidate_sites = specs.candidate_sites();

  // for each orbit generator of size n
  for (auto orbit_generator_it = begin; orbit_generator_it != end;
       ++orbit_generator_it) {
    Index orig_size = generators.elements.size();

    // print status messages
    status << clean << '\r' << "  Calculating orbit branch "
           << orbit_generator_it->size() + 1 << ":  Expanding orbit "
           << std::distance(begin, orbit_generator_it) << " / "
           << std::distance(begin, end) << "  of branch "
           << orbit_generator_it->size() << "." << std::flush;

    // by looping over each site in the grid,
    for (auto site_it = candidate_sites.first;
         site_it != candidate_sites.second; ++site_it) {
      // don't duplicate sites in cluster
      if (contains(*orbit_generator_it, *site_it)) {
        continue;
      }

      // create a test cluster from prototype
      cluster_type test(*orbit_generator_it);

      // add the new site
      test.elements().push_back(*site_it);
      // filter clusters
      if (!filter(test)) {
        continue;
      }

      // try inserting test (only uniques will be kept)
      generators.insert(test);
    }

    status << "  New orbits: " << generators.elements.size() - orig_size
           << std::endl;
  }

  return generators;
}
}  // namespace

// ------- Generating asymmetric unit orbits ---------------------------------

/// \brief Generate the asymmetric unit, using OrbitBranchSpecs
///
/// \param specs OrbitBranchSpecs for the asymmetric unit
/// \param result An output iterator for orbits of IntegralCluster
///
/// \relates IntegralCluster
///
template <typename OrbitType, typename OrbitOutputIterator>
OrbitOutputIterator make_asymmetric_unit(
    const OrbitBranchSpecs<OrbitType> &specs, OrbitOutputIterator result,
    std::ostream &status) {
  IntegralCluster empty(specs.prim());
  const SymGroup &g = specs.generating_group();

  // generate the null cluster orbit
  OrbitType null_cluster_orbit(empty, g, specs.sym_compare());

  // Use it to generate the first orbit branch
  std::vector<OrbitType> orbits(1, null_cluster_orbit);
  return make_next_orbitbranch(orbits.cbegin(), orbits.cend(), specs, result,
                               status);
}

/// \brief Use orbits of size n to generate orbits of size n+1
///
/// \param begin,end A range of input orbits of size n
/// \param specs OrbitBranchSpecs for orbits of size n+1
/// \param result An output iterator for orbits of IntegralCluster
/// \param stutus Stream for status messages
///
/// \ingroup IntegralCluster
///
template <typename OrbitType, typename OrbitInputIterator,
          typename OrbitOutputIterator>
OrbitOutputIterator make_next_orbitbranch(
    OrbitInputIterator begin, OrbitInputIterator end,
    const OrbitBranchSpecs<OrbitType> &specs, OrbitOutputIterator result,
    std::ostream &status) {
  /// Construct an OrbitGenerators object to collect orbit generating elements
  OrbitGenerators<OrbitType> generators(specs.generating_group(),
                                        specs.sym_compare());

  /// Use OrbitBranchSpecs to insert orbit generating elements for the next
  /// orbitbranch
  _insert_next_orbitbranch_generators(prototype_iterator(begin),
                                      prototype_iterator(end), specs,
                                      generators, status);

  /// Generate orbits from the orbit generating elements
  return generators.make_orbits(result);
}

/// \brief Generate Orbit<IntegralCluster> using OrbitBranchSpecs
///
/// \param begin,end A range of OrbitBranchSpecs, starting with the null branch
/// \param custom_generators A vector of custom orbit generating clusters
/// \param result An output iterator for Orbit
/// \param status Stream for status messages
///
template <typename OrbitBranchSpecsIterator, typename OrbitOutputIterator>
OrbitOutputIterator make_orbits(
    OrbitBranchSpecsIterator begin, OrbitBranchSpecsIterator end,
    const std::vector<IntegralClusterOrbitGenerator> &custom_generators,
    OrbitOutputIterator result, std::ostream &status) {
  if (begin == end) {
    throw libcasm_runtime_error(
        "Error in make_orbits: No OrbitBranchSpecs (begin==end)");
  }

  typedef typename OrbitOutputIterator::container_type container_type;
  typedef typename container_type::value_type orbit_type;

  // -- Collect orbit generating elements for each branch, to generate the next
  // branch
  std::vector<OrbitGenerators<orbit_type>> generators;
  generators.reserve(std::distance(begin, end));

  // -- Combines all orbit branches
  OrbitGenerators<orbit_type> all_generators(begin->generating_group(),
                                             begin->sym_compare());

  // branches should be ordered already, so insert with end as hint
  auto insert_branch = [&](OrbitGenerators<orbit_type> &all,
                           OrbitGenerators<orbit_type> &branch) {
    for (auto it = branch.elements.begin(); it != branch.elements.end(); ++it) {
      all.elements.insert(all.elements.end(), *it);
    }
  };
  // -- construct null cluster orbit
  status << "  Adding null orbit branch." << std::flush;

  auto specs_it = begin;
  if (specs_it != end) {
    generators.emplace_back(begin->generating_group(), begin->sym_compare());
    _insert_null_cluster_generator(specs_it->prim(), generators.back());
    ++specs_it;
    insert_branch(all_generators, generators.back());
  }
  status << "  New orbits: 1" << std::endl;

  // -- construct additional branches
  // print status messages
  std::string clean(100, ' ');

  // generate orbit branches 1+ using the previously generated branch:
  auto prev_gen = generators.begin();
  while (specs_it != end) {
    generators.emplace_back(specs_it->generating_group(),
                            specs_it->sym_compare());
    _insert_next_orbitbranch_generators(prev_gen->elements.begin(),
                                        prev_gen->elements.end(), *specs_it,
                                        generators.back(), status);
    ++specs_it;
    ++prev_gen;
    insert_branch(all_generators, generators.back());
  }

  // -- add custom orbit generators

  for (int i = 0; i < custom_generators.size(); ++i) {
    status << "  Adding custom orbit " << i << "/" << custom_generators.size()
           << "." << std::flush;
    if (custom_generators[i].include_subclusters) {
      status << " Include subclusters." << std::flush;
    }

    Index orig_size = all_generators.elements.size();

    all_generators.insert(custom_generators[i].prototype);
    if (custom_generators[i].include_subclusters) {
      insert_subcluster_generators(custom_generators[i].prototype,
                                   all_generators);
    }

    status << "  New orbits: " << all_generators.elements.size() - orig_size
           << std::endl;
  }

  status << std::endl;

  // make orbits
  return all_generators.make_orbits(result);
}

/* -- SymCompareType-Specific IntegralCluster Orbit functions --- */

/// \brief Generate the prim periodic asymmetric unit
///
/// \param prim A PrimType
/// \param site_filter A filter function that returns true for Site that
///        should be considered for the neighborhood (i.e. to check the number
///        of components)
/// \param xtal_tol Crystallography tolerance
/// \param result An output iterator for orbits of IntegralCluster
/// \param status Stream for status messages
///
/// - Uses prim.factor_group as the generating group
/// - Uses PrimPeriodicSymCompare<IntegralCluster>(xtal_tol) for cluster
/// equivalence
/// - Figures out candidate_sites from max_length and site_filter input to
///   create OrbitBranchSpecs and calls make_orbits
///
/// \relates IntegralCluster
///
template <typename OrbitOutputIterator>
OrbitOutputIterator make_prim_periodic_asymmetric_unit(
    std::shared_ptr<Structure const> prim_ptr,
    SiteFilterFunction const &site_filter, double xtal_tol,
    OrbitOutputIterator result, std::ostream &status) {
  typedef PrimPeriodicOrbit<IntegralCluster> orbit_type;
  typedef typename orbit_type::Element cluster_type;
  typedef PrimPeriodicSymCompare<IntegralCluster> symcompare_type;

  SymGroup const &generating_grp = prim_ptr->factor_group();

  symcompare_type sym_compare(prim_ptr, xtal_tol);

  std::vector<xtal::UnitCellCoord> candidate_sites;
  for (int i = 0; i < prim_ptr->basis().size(); ++i) {
    if (site_filter(prim_ptr->basis()[i])) {
      candidate_sites.emplace_back(i, 0, 0, 0);
    }
  }

  OrbitBranchSpecs<orbit_type> specs(
      *prim_ptr, candidate_sites.begin(), candidate_sites.end(), generating_grp,
      [](cluster_type const &test) { return true; }, sym_compare);

  return make_asymmetric_unit(specs, result, status);
}

/// \brief Generate Orbit<IntegralCluster> by specifying max cluster length for
/// each branch
///
/// \param prim Primitive structure
/// \param max_length vector of max_length of pairs of cluster sites. Expects
///        that max_length[b] is the max_length for orbit branch b. The values
///        for the null cluster and point clusters are ignored.
/// \param custom_generators A vector of custom orbit generating clusters
/// \param site_filter A filter function that returns true for Site that
///        should be considered for the neighborhood (i.e. to check the number
///        of components)
/// \param xtal_tol Crystallography tolerance
/// \param result An output iterator for Orbit
/// \param status Stream for status messages
///
/// - Uses prim.factor_group as the generating group
/// - Uses PrimPeriodicSymCompare<IntegralCluster>(xtal_tol) for cluster
/// equivalence
/// - Figures out candidate_sites from max_length and site_filter input to
///   create OrbitBranchSpecs and calls make_orbits
///
/// \relates IntegralCluster
template <typename OrbitOutputIterator>
OrbitOutputIterator make_prim_periodic_orbits(
    std::shared_ptr<Structure const> prim_ptr,
    std::vector<double> const &max_length,
    std::vector<IntegralClusterOrbitGenerator> const &custom_generators,
    SiteFilterFunction const &site_filter, double xtal_tol,
    OrbitOutputIterator result, std::ostream &status) {
  typedef PrimPeriodicOrbit<IntegralCluster> orbit_type;
  typedef typename orbit_type::Element cluster_type;
  typedef PrimPeriodicSymCompare<IntegralCluster> symcompare_type;

  SymGroup const &generating_grp = prim_ptr->factor_group();

  symcompare_type sym_compare(prim_ptr, xtal_tol);

  // collect OrbitBranchSpecs here
  std::vector<OrbitBranchSpecs<orbit_type>> specs;

  // collect the environment of sites here
  std::vector<xtal::UnitCellCoord> candidate_sites;

  // --- add specs for null cluster orbit ------------------
  if (max_length.size() >= 1) {
    specs.emplace_back(
        *prim_ptr, candidate_sites.begin(), candidate_sites.end(),
        generating_grp, [](cluster_type const &test) { return true; },
        sym_compare);
  }

  // --- add specs for asymmetric unit orbit ------------------
  if (max_length.size() >= 2) {
    for (int i = 0; i < prim_ptr->basis().size(); ++i) {
      if (site_filter(prim_ptr->basis()[i])) {
        candidate_sites.emplace_back(i, 0, 0, 0);
      }
    }
    specs.emplace_back(
        *prim_ptr, candidate_sites.begin(), candidate_sites.end(),
        generating_grp, [](cluster_type const &test) { return true; },
        sym_compare);
  }

  // --- add specs for additional orbit branches ------------------
  for (auto it = max_length.begin() + 2; it != max_length.end(); ++it) {
    // construct the neighborhood of sites to consider for the orbit
    candidate_sites.clear();
    neighborhood(*prim_ptr, *it, site_filter,
                 std::back_inserter(candidate_sites), xtal_tol);

    auto max_length_filter = [=](cluster_type const &test) {
      return sym_compare.make_invariants(test).displacement().back() < *it;
    };

    specs.emplace_back(*prim_ptr, candidate_sites.begin(),
                       candidate_sites.end(), generating_grp, max_length_filter,
                       sym_compare);
  }

  // now generate orbits
  return make_orbits(specs.begin(), specs.end(), custom_generators, result,
                     status);
}

/// \brief Iterate over all sites in an orbit and insert a UnitCellCoord
///
/// \param orbit an Orbit<IntegralCluster>
/// \param result an OutputIterator for UnitCellCoord
///
/// \result the resulting OutputIterator
///
/// This simply outputs all UnitCellCoord for clusters that include the origin
/// UnitCell, without any standard order. It uses all clusters that touch origin
/// unitcell, including translationally equivalent clusters.
///
/// \ingroup IntegralCluster
///
template <typename OutputIterator>
OutputIterator prim_periodic_orbit_neighborhood(
    PrimPeriodicOrbit<IntegralCluster> const &orbit, OutputIterator result) {
  for (const auto &equiv : orbit) {
    // UnitCellCoord for all sites in cluster
    std::vector<xtal::UnitCellCoord> coord(equiv.begin(), equiv.end());

    // UnitCellCoord for 'flowertree': all clusters that touch origin unitcell
    //  (includes translationally equivalent clusters)
    for (int ns_i = 0; ns_i < coord.size(); ++ns_i) {
      for (int ns_j = 0; ns_j < coord.size(); ++ns_j) {
        *result++ = xtal::UnitCellCoord(
            coord[ns_j].sublattice(),
            coord[ns_j].unitcell() - coord[ns_i].unitcell());
      }
    }
  }
  return result;
}

/// \brief Iterate over all sites in all orbits and insert a UnitCellCoord
///
/// \param begin,end Range of Orbit<IntegralCluster>
/// \param result an OutputIterator for UnitCellCoord
///
/// This simply outputs all UnitCellCoord for clusters that include the origin
/// UnitCell, without any standard order. It uses all clusters that touch origin
/// unitcell, including translationally equivalent clusters.
///
/// \ingroup IntegralCluster
///
template <typename ClusterOrbitIterator, typename OutputIterator>
OutputIterator prim_periodic_neighborhood(ClusterOrbitIterator begin,
                                          ClusterOrbitIterator end,
                                          OutputIterator result) {
  // create a neighborhood of all UnitCellCoord that an Orbitree touches
  for (auto it = begin; it != end; ++it) {
    result = prim_periodic_orbit_neighborhood(*it, result);
  }
  return result;
}

/// \brief Iterate over all sites in an orbit and insert a UnitCellCoord
///
/// \param orbit an PrimPeriodicOrbit<IntegralCluster>
/// \param result an OutputIterator for UnitCellCoord
///
/// \result the resulting OutputIterator
///
/// This simply outputs all UnitCellCoord for clusters that include the origin
/// UnitCell, without any standard order. It uses all clusters that touch origin
/// unitcell, including translationally equivalent clusters.
///
/// \ingroup IntegralCluster
///
template <typename OutputIterator, typename OrbitType>
OutputIterator flower_neighborhood(OrbitType const &orbit,
                                   OutputIterator result) {
  xtal::UnitCellCoord const *ucc_ptr(nullptr);
  for (auto const &equiv : orbit) {
    for (xtal::UnitCellCoord const &ucc : equiv) {
      ucc_ptr = &ucc;
      break;
    }
    if (ucc_ptr) break;
  }

  if (!ucc_ptr) return result;

  SymGroup identity_group(
      orbit.prototype().prim().factor_group().begin(),
      (orbit.prototype().prim().factor_group().begin()) + 1);
  OrbitType empty_orbit(typename OrbitType::Element(orbit.prototype().prim()),
                        identity_group, orbit.sym_compare());
  typename OrbitType::Element test(empty_orbit.prototype());
  test.elements().push_back(*ucc_ptr);

  // Loop over each site in each cluster of the orbit
  for (auto const &equiv : orbit) {
    for (xtal::UnitCellCoord const &ucc : equiv) {
      // create a test cluster from prototype
      // add the new site
      test.elements()[0] = ucc;

      test = orbit.sym_compare().prepare(test);

      xtal::UnitCell trans = test.element(0).unitcell() - ucc.unitcell();
      for (xtal::UnitCellCoord const &ucc2 : equiv) {
        *result++ = (ucc2 + trans);
      }
    }
  }

  return result;
}

/// \brief Iterate over all sites in all orbits and insert a UnitCellCoord
///
/// \param begin,end Range of PrimPeriodicOrbit<IntegralCluster>
/// \param result an OutputIterator for UnitCellCoord
///
/// This simply outputs all UnitCellCoord for clusters that include the origin
/// UnitCell, without any standard order. It uses all clusters that touch origin
/// unitcell, including translationally equivalent clusters.
///
/// \ingroup IntegralCluster
///
template <typename ClusterOrbitIterator, typename OutputIterator>
OutputIterator flower_neighborhood(ClusterOrbitIterator begin,
                                   ClusterOrbitIterator end,
                                   OutputIterator result) {
  // create a neighborhood of all UnitCellCoord that an Orbitree touches
  for (auto it = begin; it != end; ++it) {
    result = flower_neighborhood(*it, result);
  }
  return result;
}

/// \brief Return index of asymmetric unit containing unitcellcoord
///
/// \param begin,end Range of orbits of IntegralCluster
/// \param unitcellcoord Site to search for
///
/// \returns Index of point orbit (starting from 0, counting point orbits only)
/// that contains the point cluster consisting of unitcellcoord. Uses
/// sym_compare to prepare the point cluster. Returns -1 if not found in any
/// point orbit.
///
template <typename ClusterOrbitIterator>
Index find_asymmetric_unit_index(xtal::UnitCellCoord const &unitcellcoord,
                                 ClusterOrbitIterator begin,
                                 ClusterOrbitIterator end) {
  if (begin == end) {
    return -1;
  }

  auto orbit_it = begin;
  auto const &prim = orbit_it->prototype().prim();
  auto const &sym_compare = orbit_it->sym_compare();

  // create a test cluster with just unitcellcoord, then "prepare" it for
  // comparisons
  IntegralCluster test{prim};
  test.elements().push_back(unitcellcoord);
  test = sym_compare.prepare(test);

  // find point orbit that contains unitcellcoord
  Index asym_unit_index = 0;
  for (; orbit_it != end; ++orbit_it) {
    if (orbit_it->prototype().size() != 1) {
      continue;
    }
    auto cluster_it = orbit_it->find(test);
    if (cluster_it != orbit_it->end()) {
      return asym_unit_index;
    }
    ++asym_unit_index;
  }

  // if unitcellcoord not included in any point orbit
  return -1;
}

}  // namespace CASM

#endif
