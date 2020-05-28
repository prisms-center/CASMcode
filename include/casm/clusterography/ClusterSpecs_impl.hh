#ifndef CASM_ClusterSpecs_impl
#define CASM_ClusterSpecs_impl

#include "casm/clusterography/ClusterSpecs.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clusterography/IntegralCluster_impl.hh"
#include "casm/crystallography/Site.hh"
#include "casm/symmetry/Orbit.hh"

namespace CASM {

  /* -- Generate prim periodic orbits --------------------------------------- */

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
  /// - Uses PrimPeriodicSymCompare<IntegralCluster>(xtal_tol) for cluster equivalence
  /// - Figures out candidate_sites from max_length and site_filter input to
  ///   create OrbitBranchSpecs and calls make_orbits
  ///
  /// \relates IntegralCluster
  ///
  template<typename OrbitOutputIterator>
  OrbitOutputIterator make_prim_periodic_asymmetric_unit(
    std::shared_ptr<Structure const> prim_ptr,
    SiteFilterFunction const &site_filter,
    double xtal_tol,
    OrbitOutputIterator result,
    std::ostream &status) {

    typedef PrimPeriodicOrbit<IntegralCluster> orbit_type;
    typedef typename orbit_type::Element cluster_type;
    typedef PrimPeriodicSymCompare<IntegralCluster> symcompare_type;

    SymGroup const &generating_grp = prim_ptr->factor_group();

    symcompare_type sym_compare(prim_ptr, xtal_tol);

    std::vector<xtal::UnitCellCoord> candidate_sites;
    for(int i = 0; i < prim_ptr->basis().size(); ++i) {
      if(site_filter(prim_ptr->basis()[i])) {
        candidate_sites.emplace_back(i, 0, 0, 0);
      }
    }

    OrbitBranchSpecs<orbit_type> specs(
      *prim_ptr,
      candidate_sites.begin(),
      candidate_sites.end(),
      generating_grp,
    [](cluster_type const & test) {
      return true;
    },
    sym_compare);

    return make_asymmetric_unit(specs, result, status);
  }

  /// \brief Generate Orbit<IntegralCluster> by specifying max cluster length for each branch
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
  /// - Uses PrimPeriodicSymCompare<IntegralCluster>(xtal_tol) for cluster equivalence
  /// - Figures out candidate_sites from max_length and site_filter input to
  ///   create OrbitBranchSpecs and calls make_orbits
  ///
  /// \relates IntegralCluster
  template<typename OrbitOutputIterator>
  OrbitOutputIterator make_prim_periodic_orbits(
    std::shared_ptr<Structure const> prim_ptr,
    std::vector<double> const &max_length,
    std::vector<IntegralClusterOrbitGenerator> const &custom_generators,
    SiteFilterFunction const &site_filter,
    double xtal_tol,
    OrbitOutputIterator result,
    std::ostream &status) {
    typedef PrimPeriodicOrbit<IntegralCluster> orbit_type;
    typedef typename orbit_type::Element cluster_type;
    typedef PrimPeriodicSymCompare<IntegralCluster> symcompare_type;

    SymGroup const &generating_grp = prim_ptr->factor_group();

    symcompare_type sym_compare(prim_ptr, xtal_tol);

    // collect OrbitBranchSpecs here
    std::vector<OrbitBranchSpecs<orbit_type> > specs;

    // collect the environment of sites here
    std::vector<xtal::UnitCellCoord> candidate_sites;

    // --- add specs for null cluster orbit ------------------
    if(max_length.size() >= 1) {
      specs.emplace_back(*prim_ptr,
                         candidate_sites.begin(),
                         candidate_sites.end(),
                         generating_grp,
      [](cluster_type const & test) {
        return true;
      },
      sym_compare);
    }


    // --- add specs for asymmetric unit orbit ------------------
    if(max_length.size() >= 2) {
      for(int i = 0; i < prim_ptr->basis().size(); ++i) {
        if(site_filter(prim_ptr->basis()[i])) {
          candidate_sites.emplace_back(i, 0, 0, 0);
        }
      }
      specs.emplace_back(*prim_ptr,
                         candidate_sites.begin(),
                         candidate_sites.end(),
                         generating_grp,
      [](cluster_type const & test) {
        return true;
      },
      sym_compare);
    }

    // --- add specs for additional orbit branches ------------------
    for(auto it = max_length.begin() + 2; it != max_length.end(); ++it) {

      // construct the neighborhood of sites to consider for the orbit
      candidate_sites.clear();
      neighborhood(*prim_ptr, *it, site_filter, std::back_inserter(candidate_sites), xtal_tol);

      auto max_length_filter = [ = ](cluster_type const & test) {
        return sym_compare.make_invariants(test).displacement().back() < *it;
      };

      specs.emplace_back(*prim_ptr,
                         candidate_sites.begin(),
                         candidate_sites.end(),
                         generating_grp,
                         max_length_filter,
                         sym_compare);
    }

    // now generate orbits
    return make_orbits(specs.begin(), specs.end(), custom_generators, result, status);

  }
}

#endif
