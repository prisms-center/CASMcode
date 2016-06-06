#ifndef CASM_ClusterOrbits_impl
#define CASM_ClusterOrbits_impl

#include "casm/clusterography/ClusterOrbits.hh"

namespace CASM {

  /* -- Cluster Orbit generating function definitions ------------------------------------- */


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

    std::vector<bool> assigned(prim.basis.size(), false);

    // for all sites in the basis
    for(int i = 0; i < assigned.size(); i++) {
      if(assigned[i])
        continue;

      // create a prototype cluster
      IntegralCluster cluster(prim);

      cluster.sites().push_back(UnitCellCoord(prim, UnitCell(0, 0, 0), i));

      // generate an orbit
      OrbitType orbit(cluster, generating_grp, sym_compare);

      // note sites that have been added to the orbit
      for(int j = 0; j < orbit.size(); j++) {
        assigned[orbit[j][0].sublat()] = true;
      }

      // add the orbit to the asym_unit
      *result++ = orbit;
    }

    return result;
    `

  }

  /// \brief Generate the asymmetric unit, using OrbitBranchSpecs
  ///
  /// \param specs OrbitBranchSpecs for the asymmetric unit
  /// \param result An output iterator for orbits of IntegralCluster
  ///
  /// \relates IntegralCluster
  ///
  template<typename OrbitType, typename OrbitOutputIterator>
  OrbitOutputIterator make_asymmetric_unit(const OrbitBrancSpecs<OrbitType> &specs, OrbitOutputIterator result) {

    IntegralCluster empty(specs.prim());
    const SymGroup &g = specs.generating_group();

    // generate the null cluster orbit
    null_cluster_orbit = OrbitType(empty, g.begin(), g.end(), specs.sym_compare());

    // Use it to generate the first orbit branch
    std::vector<OrbitType> orbits(1, null_cluster_orbit);
    return next_orbitbranch(orbits.cbegin(), orbits.cend(), specs, result);
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
    //std::cout << "begin next_orbitbranch" << std::endl;

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
    for(const auto &val : orbits) {
      result = std::move(val.second.begin(), val.second.end(), result);
    }

    // print status messages
    status << clean << '\r';

    //std::cout << "finish next_orbitbranch" << std::endl;
    return result;
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

    typedef typename std::iterator_traits<OrbitBranchSpecsIterator>::value_type::OrbitType OrbitType;

    // Temporarily store Orbits because we need to use ranges of them to
    //   construct successive orbit branches
    std::vector<OrbitType> _orbits;

    // -- construct null cluster orbit

    //std::cout << "begin orbits()" << std::endl;
    auto specs_it = begin;

    // generate the null cluster orbit
    IntegralCluster empty(specs_it->prim());
    null_cluster_orbit = OrbitType(empty, specs_it->generating_group(), specs_it->sym_compare());


    // -- construct asymmetric unit

    // print status message
    std::string clean(80, ' ');
    status << clean << '\r' << "Calculating asymmetric unit\r" << std::flush;

    ++specs_it;

    // generate asymmetic unit, store iterators to begin and end of size=1 orbits
    asymmetric_unit(*specs_it, std::back_inserter(_orbits));


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
      next_orbitbranch(_orbits.cbegin() + a, _orbits.cend(), *specs_it, std::back_inserter(_orbits));

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

    // add specs for null cluster orbit
    specs.emplace_back(prim,
                       candidate_sites.begin(),
                       candidate_sites.end()
                       generating_grp,
    [](const cluster_type & test) {
      return true;
    },
    sym_compare);


    // add specs for asymmetric unit orbit
    for(int i = 0; i < prim.basis.size(); ++i) {
      candidate_sites.emplace_back(prim, i, 0, 0, 0);
    }
    specs.emplace_back(prim,
                       candidate_sites.begin(),
                       candidate_sites.end()
                       generating_grp,
    [](const cluster_type & test) {
      return true;
    },
    sym_compare);

    // add specs for additional orbit branches
    for(auto it = max_length.begin(), it != max_length.end(); ++it) {

      // construct the neighborhood of sites to consider for the orbit
      candidate_sites.clear();
      neighborhood(prim, *it, site_filter, std::back_inserter(candidate_sites), xtal_tol);

      specs.emplace_back(prim,
                         candidate_sites.begin(),
                         candidate_sites.end()
                         generating_grp,
      [](const cluster_type & test) {
        return test.invariants().displacement().back() < *it;
      },
      sym_compare);
    }

    // now generate orbits
    return orbits(specs.begin(), specs.end(), result, status);
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
    const std::function<bool (Site)> &site_filter,
    const SymCompareType &sym_compare,
    OrbitOutputIterator result,
    std::ostream &status) {

    typedef Orbit<IntegralCluster, SymCompareType> orbit_type;

    // get 'max_length' vector from bspecs
    std::vector<double> max_length;

    auto update_max_length = [&](int branch, double length) {
      while(branch > max_length.size() - 1) {
        max_length.push_back(0.0);
      }
      max_length[branch] = length;
    };

    const auto &j = json["orbit_branch_specs"];
    for(auto it = j.begin(); it != j.end(); ++it) {
      update_max_length(std::stoi(it.name()), it->find("max_length")->get<double>());
    }

    // generate clusters
    std::vector<orbit_type> orbits;
    make_orbits(prim, generating_grp, max_length, site_filter, sym_compare, std::back_inserter(orbits), status);

    // read custom orbit specs
    // generate custom clusters
    // make_custom_orbits(...)
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
                                                       unsigned int size) {

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