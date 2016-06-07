#ifndef CASM_ClusterOrbits_impl
#define CASM_ClusterOrbits_impl

#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/misc/algorithm.hh"

namespace CASM {

  namespace {

    /// Read max_length vector from 'bspecs' JSON
    ///
    /// \returns std::vector<double> giving 'max_length' for clusters in branch 2+ (pairs, triplets, etc.)
    ///
    /// probably should get organized somewhere...
    std::vector<double> max_length_from_bspecs(const jsonParser &bspecs) {

      std::vector<double> max_length;

      auto update_max_length = [&](int branch, double length) {
        while(branch - 1 > max_length.size()) {
          max_length.push_back(0.0);
        }
        max_length[branch - 2] = length;
      };

      if(bspecs.contains("orbit_branch_specs")) {
        const auto &j = bspecs["orbit_branch_specs"];
        for(auto it = j.begin(); it != j.end(); ++it) {
          update_max_length(std::stoi(it.name()), it->find("max_length")->get<double>());
        }
      }

      return max_length;
    }

  }

  /* -- Cluster Orbit generating function definitions ------------------------------------- */

  /// \brief Output the neighborhood of UnitCellCoord within max_radius of a unit cell
  ///
  /// \param unit The unit cell Structure
  /// \param max_radius The neighborhood distance cutoff
  /// \param site_filter A filter function that returns true for CoordType that
  ///        should be considered for the neighborhood
  /// \param result Output iterator for container of UnitCellCoord
  /// \param xtal_tol Crystallography tolerance used to contstruct UnitCellCoord from CoordType
  ///
  /// \returns Output iterator after generating the neighborhood
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename CoordType, typename OutputIterator>
  OutputIterator neighborhood(
    const Structure &unit,
    double max_radius,
    std::function<bool (CoordType)> site_filter,
    OutputIterator result,
    double xtal_tol) {

    auto dim = unit.lattice().enclose_sphere(max_radius);
    EigenCounter<Eigen::Vector3i> grid_count(-dim, dim, Eigen::Vector3i::Constant(1));
    Coordinate lat_point(unit.lattice());
    const auto &basis = unit.basis;

    do {
      lat_point.frac() = grid_count().cast<double>();

      for(auto it = basis.begin(); it != basis.end(); ++it) {
        if(!site_filter(*it)) {
          continue;
        }

        Coordinate test(*it + lat_point);
        if(std::any_of(basis.begin(),
                       basis.end(),
        [&](const Coordinate & coord) {
        return test.dist(coord) < max_radius;
        })) {
          *result++ = UnitCellCoord(unit, test, xtal_tol);
        }
      }
    }
    while(++grid_count);
    return result;
  }

  /// \brief Generate the asymmetric unit, including all sites
  ///
  /// \param prim A PrimType
  /// \param generating_grp Iterators over SymOp to use for generating the asymmetric unit orbit
  /// \param compare Determines the symmetry properties used to generated the orbits
  /// \param result An output iterator for orbits of IntegralCluster
  ///
  /// \relates IntegralCluster
  ///
  template<typename OrbitOutputIterator, typename SymCompareType>
  OrbitOutputIterator make_asymmetric_unit(
    const IntegralCluster::PrimType &prim,
    const SymGroup &generating_grp,
    const SymCompareType &sym_compare,
    OrbitOutputIterator result) {

    typedef typename OrbitOutputIterator::container_type container_type;
    typedef typename container_type::value_type orbit_type;

    std::vector<bool> assigned(prim.basis.size(), false);

    // for all sites in the basis
    for(int i = 0; i < assigned.size(); i++) {
      if(assigned[i])
        continue;

      // create a prototype cluster
      IntegralCluster cluster(prim);

      cluster.elements().push_back(UnitCellCoord(prim, i, UnitCell(0, 0, 0)));

      // generate an orbit
      orbit_type orbit(cluster, generating_grp, sym_compare);

      // note sites that have been added to the orbit
      for(int j = 0; j < orbit.size(); j++) {
        assigned[orbit[j][0].sublat()] = true;
      }

      // add the orbit to the asym_unit
      *result++ = orbit;
    }

    return result;
  }

  /// \brief Generate the asymmetric unit, using OrbitBranchSpecs
  ///
  /// \param specs OrbitBranchSpecs for the asymmetric unit
  /// \param result An output iterator for orbits of IntegralCluster
  ///
  /// \relates IntegralCluster
  ///
  template<typename OrbitType, typename OrbitOutputIterator>
  OrbitOutputIterator make_asymmetric_unit(
    const OrbitBranchSpecs<OrbitType> &specs,
    OrbitOutputIterator result,
    std::ostream &status) {

    IntegralCluster empty(specs.prim());
    const SymGroup &g = specs.generating_group();

    // generate the null cluster orbit
    OrbitType null_cluster_orbit(empty, g, specs.sym_compare());

    // Use it to generate the first orbit branch
    std::vector<OrbitType> orbits(1, null_cluster_orbit);
    return make_next_orbitbranch(orbits.cbegin(), orbits.cend(), specs, result, status);
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
  template<typename OrbitType, typename OrbitInputIterator, typename OrbitOutputIterator>
  OrbitOutputIterator make_next_orbitbranch_impl0(OrbitInputIterator begin,
                                                  OrbitInputIterator end,
                                                  const OrbitBranchSpecs<OrbitType> &specs,
                                                  OrbitOutputIterator result,
                                                  std::ostream &status) {

    typedef IntegralCluster cluster_type;
    typedef cluster_type::InvariantsType invariants_type;

    const auto &sym_compare = specs.sym_compare();
    const auto &filter = specs.filter();
    const auto &g = specs.generating_group();

    // store orbits as we find them
    std::set<OrbitType> orbits;

    // print status messages
    std::string clean(100, ' ');

    // contains a pair of iterators over candidate UnitCellCoord
    auto candidate_sites = specs.candidate_sites();

    // for each orbit of size n
    for(auto orbit_it = begin; orbit_it != end; ++orbit_it) {

      // print status messages
      status << clean << '\r'
             << "  Calculating orbit branch " << orbit_it->prototype().size() + 1
             << ":  Expanding orbit " << std::distance(begin, orbit_it)
             << " / " << std::distance(begin, end)
             << "  of branch " << orbit_it->prototype().size()
             << ".  New orbits: " << orbits.size() << "\r" << std::flush;

      // by looping over each site in the grid,
      for(auto site_it = candidate_sites.first; site_it != candidate_sites.second; ++site_it) {

        // don't duplicate sites in cluster
        if(contains(orbit_it->prototype(), *site_it)) {
          continue;
        }

        // create a test cluster from prototype
        cluster_type test(orbit_it->prototype());

        // add the new site
        test.elements().push_back(*site_it);

        // 'prepare' the test cluster for comparison
        test = sym_compare.prepare(test);

        // filter clusters
        if(!filter(test)) {
          continue;
        }

        // try to find test cluster in already found orbits
        auto it = find_orbit(orbits.begin(), orbits.end(), test);
        if(it != orbits.end()) {
          continue;
        }

        // if not yet found, use test to generate a new Orbit
        orbits.insert(OrbitType(test, g, sym_compare));

      }
    }

    // output Orbits
    result = std::move(orbits.begin(), orbits.end(), result);

    // print status messages
    status << clean << '\r';

    return result;
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
  template<typename OrbitType, typename OrbitInputIterator, typename OrbitOutputIterator>
  OrbitOutputIterator make_next_orbitbranch_impl1(OrbitInputIterator begin,
                                                  OrbitInputIterator end,
                                                  const OrbitBranchSpecs<OrbitType> &specs,
                                                  OrbitOutputIterator result,
                                                  std::ostream &status) {

    typedef IntegralCluster cluster_type;
    typedef OrbitType orbit_type;
    typedef cluster_type::InvariantsType invariants_type;

    const auto &sym_compare = specs.sym_compare();
    const auto &filter = specs.filter();
    const auto &g = specs.generating_group();

    auto compare = [&](const cluster_type & A, const cluster_type & B) {
      if(A.size() != B.size()) {
        return A.size() < B.size();
      }
      return sym_compare.intra_orbit_compare(A, B);
    };

    auto canonical_equiv = [&](const cluster_type & clust) {
      cluster_type result = sym_compare.prepare(clust);
      for(const auto &op : g) {
        auto test = sym_compare.prepare(copy_apply(op, clust));
        if(sym_compare.intra_orbit_compare(result, test)) {
          result = test;
        }
      }
      return result;
    };

    // store generating elements as we find them
    std::set<IntegralCluster, decltype(compare)> generators(compare);

    // print status messages
    std::string clean(100, ' ');

    // contains a pair of iterators over candidate UnitCellCoord
    auto candidate_sites = specs.candidate_sites();

    // for each orbit of size n
    for(auto orbit_it = begin; orbit_it != end; ++orbit_it) {

      // print status messages
      status << clean << '\r'
             << "  Calculating orbit branch " << orbit_it->prototype().size() + 1
             << ":  Expanding orbit " << std::distance(begin, orbit_it)
             << " / " << std::distance(begin, end)
             << "  of branch " << orbit_it->prototype().size()
             << ".  New orbits: " << generators.size() << "\r" << std::flush;

      // by looping over each site in the grid,
      for(auto site_it = candidate_sites.first; site_it != candidate_sites.second; ++site_it) {

        // don't duplicate sites in cluster
        if(contains(orbit_it->prototype(), *site_it)) {
          continue;
        }

        // create a test cluster from prototype
        cluster_type test(orbit_it->prototype());

        // add the new site
        test.elements().push_back(*site_it);

        // 'prepare' the test cluster for comparison
        test = sym_compare.prepare(test);

        // put into a 'canonical' equivalent form
        test = canonical_equiv(test);

        // filter clusters
        if(!filter(test)) {
          continue;
        }

        // try inserting test (only uniques will be kept)
        generators.insert(test);
      }
    }

    // generate sorted orbits
    std::set<orbit_type> orbits;
    for(const auto &e : generators) {
      orbits.insert(orbit_type(e, g, sym_compare));
    }

    // output Orbits
    result = std::move(orbits.begin(), orbits.end(), result);

    // print status messages
    status << clean << '\r';

    return result;
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
  template<typename OrbitType, typename OrbitInputIterator, typename OrbitOutputIterator>
  OrbitOutputIterator make_next_orbitbranch(OrbitInputIterator begin,
                                            OrbitInputIterator end,
                                            const OrbitBranchSpecs<OrbitType> &specs,
                                            OrbitOutputIterator result,
                                            std::ostream &status) {
    return make_next_orbitbranch_impl1(
             begin,
             end,
             specs,
             result,
             status);
  }

  /// \brief Generate Orbit using OrbitBranchSpecs
  ///
  /// \param begin,end OrbitBranchSpecsIterators over range of OrbitBranchSpecs,
  ///                  with one for each orbit branch to be calculated
  /// \param result OutputIterator for resulting Orbit
  /// \param stutus Stream for status messages
  ///
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename OrbitBranchSpecsIterator, typename OrbitOutputIterator>
  OrbitOutputIterator make_orbits(
    OrbitBranchSpecsIterator begin,
    OrbitBranchSpecsIterator end,
    OrbitOutputIterator result,
    std::ostream &status) {

    typedef typename OrbitOutputIterator::container_type container_type;
    typedef typename container_type::value_type orbit_type;

    // Temporarily store Orbits because we need to use ranges of them to
    //   construct successive orbit branches
    std::vector<orbit_type> _orbits;

    // -- construct null cluster orbit

    //std::cout << "begin orbits()" << std::endl;
    auto specs_it = begin;

    // generate the null cluster orbit
    IntegralCluster empty(specs_it->prim());
    orbit_type null_cluster_orbit(empty, specs_it->generating_group(), specs_it->sym_compare());


    // -- construct asymmetric unit

    // print status message
    std::string clean(80, ' ');
    status << clean << '\r' << "Calculating asymmetric unit\r" << std::flush;

    ++specs_it;

    // generate asymmetic unit, store iterators to begin and end of size=1 orbits
    make_asymmetric_unit(*specs_it, std::back_inserter(_orbits), status);


    // -- construct additional branches

    // generate the second orbit branch from the orbits in range [a, b), where b is the current end
    int a = 1;
    int b = _orbits.size();

    // generate second+ orbit branches:
    ++specs_it;
    for(; specs_it != end; ++specs_it) {

      // print status message
      status << clean << '\r' << "Calculating orbit branch "
             << std::distance(begin, specs_it) << "\r" << std::flush;

      // generate an orbitbranch from the orbits in range [a, next_a), where next_a is the current end
      make_next_orbitbranch(_orbits.cbegin() + a, _orbits.cend(), *specs_it, std::back_inserter(_orbits), status);

      // update the range of orbits for the next orbit branch
      a = b;
      b = _orbits.size();

    }

    // copy results
    std::copy(_orbits.cbegin(), _orbits.cend(), result);

    //std::cout << "finish orbits()" << std::endl;
    return result;
  }



  /// \brief Generate Orbit by specifying max cluster length for each branch
  ///
  /// \param prim Primitive structure
  /// \param generating_grp Iterators over SymOp to use for generating the asymmetric unit orbit
  /// \param max_length Max cluster site displacement, for branches 2+
  /// \param site_filter A filter function that returns true for Site that
  ///        should be considered for the neighborhood (i.e. to check the number of components)
  /// \param sym_compare Determines the symmetry properties used to generated the orbits
  /// \param result An output iterator for Orbit
  /// \param stutus Stream for status messages
  ///
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename OrbitOutputIterator, typename SymCompareType>
  OrbitOutputIterator make_orbits(
    const IntegralCluster::PrimType &prim,
    const SymGroup &generating_grp,
    const std::vector<double> &max_length,
    double crystallography_tol,
    const std::function<bool (Site)> &site_filter,
    const SymCompareType &sym_compare,
    OrbitOutputIterator result,
    std::ostream &status) {

    typedef IntegralCluster cluster_type;
    typedef Orbit<IntegralCluster, SymCompareType> orbit_type;

    // collect OrbitBranchSpecs here
    std::vector<OrbitBranchSpecs<orbit_type> > specs;

    // collect the environment of sites here
    std::vector<UnitCellCoord> candidate_sites;

    // --- add specs for null cluster orbit ------------------
    specs.emplace_back(prim,
                       candidate_sites.begin(),
                       candidate_sites.end(),
                       generating_grp,
    [](const cluster_type & test) {
      return true;
    },
    sym_compare);


    // --- add specs for asymmetric unit orbit ------------------
    for(int i = 0; i < prim.basis.size(); ++i) {
      candidate_sites.emplace_back(prim, i, 0, 0, 0);
    }
    specs.emplace_back(prim,
                       candidate_sites.begin(),
                       candidate_sites.end(),
                       generating_grp,
    [](const cluster_type & test) {
      return true;
    },
    sym_compare);


    // --- add specs for additional orbit branches ------------------
    for(auto it = max_length.begin(); it != max_length.end(); ++it) {

      // construct the neighborhood of sites to consider for the orbit
      candidate_sites.clear();
      neighborhood(prim, *it, site_filter, std::back_inserter(candidate_sites), crystallography_tol);

      specs.emplace_back(prim,
                         candidate_sites.begin(),
                         candidate_sites.end(),
                         generating_grp,
      [ = ](const cluster_type & test) {
        return test.invariants().displacement().back() < *it;
      },
      sym_compare);
    }

    // now generate orbits
    return make_orbits(specs.begin(), specs.end(), result, status);
  }

  /// \brief Generate Orbit from bspecs.json-type JSON input file
  ///
  /// \param prim Primitive structure
  /// \param generating_grp Iterators over SymOp to use for generating the asymmetric unit orbit
  /// \param bspecs 'bspecs.json'-like JSON object
  /// \param site_filter A filter function that returns true for Site that
  ///        should be considered for the neighborhood (i.e. to check the number of components)
  /// \param sym_compare Determines the symmetry properties used to generated the orbits
  /// \param result An output iterator for Orbit
  /// \param stutus Stream for status messages
  ///
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename OrbitOutputIterator, typename SymCompareType>
  OrbitOutputIterator make_orbits(
    const IntegralCluster::PrimType &prim,
    const SymGroup &generating_grp,
    const jsonParser &bspecs,
    double crystallography_tol,
    const std::function<bool (Site)> &site_filter,
    const SymCompareType &sym_compare,
    OrbitOutputIterator result,
    std::ostream &status) {

    typedef Orbit<IntegralCluster, SymCompareType> orbit_type;

    // generate orbits
    std::vector<orbit_type> orbits;
    make_orbits(
      prim,
      generating_grp,
      max_length_from_bspecs(bspecs),
      crystallography_tol,
      site_filter,
      sym_compare,
      std::back_inserter(orbits),
      status
    );

    // read custom orbit specs
    // generate custom clusters
    // make_custom_orbits(...)
    if(bspecs.contains("orbit_specs")) {
      throw std::runtime_error("Error: the orbit_specs generation is being re-implemented");
    }

    // output orbits
    return std::move(orbits.begin(), orbits.end(), result);
  }

  /*
    /// \brief Output all subclusters (whether symmetrically equivalent or not)
    template<typename ClusterType, typename OutputIterator>
    OutputIterator subclusters(
      const ClusterType& cluster,
      OutputIterator result) {

      // --- enumerate subclusters ----

      // count over subclusters to create test subcluster

        // prepare, and intra-orbit compare

      // ... insert unique
      // ... output
    }

    OutputIterator subcluster_orbits(
      const ClusterType& cluster,
      OutputIterator result,
      const SymGroup& grp,
      const SymCompare<ClusterType>& sym_compare);
  */

  /*
    "orbit_specs" : [
      {
        "coordinate_mode" : "Direct",
        "prototype" : [
          [ 0.000000000000, 0.000000000000, 0.000000000000 ],
          [ 1.000000000000, 0.000000000000, 0.000000000000 ],
          [ 2.000000000000, 0.000000000000, 0.000000000000 ],
          [ 3.000000000000, 0.000000000000, 0.000000000000 ]
        ],
        "include_subclusters" : true
      },
      {
        "coordinate_mode" : "Integral",
        "prototype" : [
          [ 0, 0, 0, 0 ],
          [ 1, 0, 0, 0 ],
          [ 0, 0, 0, 3 ]
        ],
        "include_subclusters" : true
      }
    ]
  */


  /* -- Cluster Orbit access/usage function definitions ------------------------------------- */


  /// \brief Returns the first range containing orbits of the requested orbit branch in the given range of Orbit
  ///
  /// \ingroup Clusterography
  ///
  template<typename OrbitIterator>
  std::pair<OrbitIterator, OrbitIterator> orbit_branch(OrbitIterator begin,
                                                       OrbitIterator end,
                                                       Index size) {

    auto branch_begin = std::find_if(begin,
                                     end,
    [&](const typename OrbitIterator::value_type & orbit) {
      return orbit.prototype().size() == size;
    });

    auto branch_end = std::find_if(branch_begin,
                                   end,
    [&](const typename OrbitIterator::value_type & orbit) {
      return orbit.prototype().size() == size + 1;
    });

    return std::pair<OrbitIterator, OrbitIterator>(branch_begin, branch_end);
  }

}

#endif