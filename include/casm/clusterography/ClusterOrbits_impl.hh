#ifndef CASM_ClusterOrbits_impl
#define CASM_ClusterOrbits_impl

#include <boost/iterator/transform_iterator.hpp>
#include <boost/lexical_cast.hpp>
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/container/Counter.hh"
#include "casm/misc/algorithm.hh"
#include "casm/misc/TypeInfo.hh"
#include "casm/misc/unique_map.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/EnumIO.hh"
#include "casm/casm_io/Log.hh"
#include "casm/symmetry/Orbit.hh"
#include "casm/symmetry/OrbitGeneration.hh"
#include "casm/symmetry/InvariantSubgroup.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clusterography/SubClusterGenerator.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clusterography/ClusterSymCompare.hh"
#include "casm/clex/Supercell.hh"
#include "casm/kinetics/DiffusionTransformation.hh"

namespace CASM {

  namespace {

    /// Read max_length vector from 'bspecs' JSON
    ///
    /// \returns std::vector<double> giving 'max_length' for clusters in branch 0, 1, etc.
    ///
    std::vector<double> max_length_from_bspecs(const jsonParser &bspecs) {

      std::vector<double> max_length;

      auto update_max_length = [&](int branch, double length) {
        while(branch >= max_length.size()) {
          max_length.push_back(0.0);
        }
        max_length[branch] = length;
      };

      if(bspecs.contains("orbit_branch_specs")) {
        const auto &j = bspecs["orbit_branch_specs"];
        for(auto it = j.begin(); it != j.end(); ++it) {
          if(it->find("max_length") != it->end()) {
            update_max_length(std::stoi(it.name()), it->find("max_length")->get<double>());
          }
        }
      }

      return max_length;
    }

    /// Read cutoff_radius vector from 'bspecs' JSON
    ///
    /// \returns std::vector<double> giving 'cutoff_radius' for clusters in branch 0, 1, etc.
    ///
    std::vector<double> cutoff_radius_from_bspecs(const jsonParser &bspecs) {

      std::vector<double> cutoff_radius;

      auto update_cutoff_radius = [&](int branch, double length) {
        while(branch >= cutoff_radius.size()) {
          cutoff_radius.push_back(0.0);
        }
        cutoff_radius[branch] = length;
      };

      if(bspecs.contains("orbit_branch_specs")) {
        const auto &j = bspecs["orbit_branch_specs"];
        for(auto it = j.begin(); it != j.end(); ++it) {
          update_cutoff_radius(std::stoi(it.name()), it->find("cutoff_radius")->get<double>());
        }
      }

      return cutoff_radius;
    }
  }

  /*
  {
    "cluster_specs": {
      "method": "ClustersByMaxLength",
      "kwargs": {
        "orbit_branch_specs": {
          "2": {
            "max_length": 6.0
          },
          "3": {
            "max_length": 6.0
          }
        },
        "orbit_specs": [
          {
            "coordinate_mode" : "Direct",
            "prototype" : [
              [ 0.000000000000, 0.000000000000, 0.000000000000 ],
              [ 1.000000000000, 0.000000000000, 0.000000000000 ],
              [ 2.000000000000, 0.000000000000, 0.000000000000 ],
              [ 3.000000000000, 0.000000000000, 0.000000000000 ]],
            "include_subclusters" : true
          },
          ...
        ]
      }
    }
  }

  {
    "cluster_specs": {
    "method": "ClustersByNumber",
    "kwargs": {
      "orbit_branch_specs": {
        "2": 20,
        "3": 100
      },
      "orbit_specs": [
        {
          "coordinate_mode" : "Direct",
          "prototype" : [
            [ 0.000000000000, 0.000000000000, 0.000000000000 ],
            [ 1.000000000000, 0.000000000000, 0.000000000000 ],
            [ 2.000000000000, 0.000000000000, 0.000000000000 ],
            [ 3.000000000000, 0.000000000000, 0.000000000000 ]],
          "include_subclusters" : true
        },
        ...]
      }
    }
  }
  */

  /*

      /// Help parsing JSON input by collecting warnings and errors
      ///
      /// Assumes input of form:
      /// \code
      /// {
      ///   "method": "MethodName",
      ///   "kwargs": {
      ///     "option1": ...,
      ///     "option2": {
      ///       "method": "SubMethodName",
      ///       "kwargs": {
      ///         "option1": ...,
      ///         "option2": ...
      ///       }
      ///     }
      ///   }
      /// }
      /// \endcode
      ///
      /// or (with top-level method already known):
      /// \code
      /// {
      ///   "option1": ...,
      ///   "option2": {
      ///     "method": "SubMethodName",
      ///     "kwargs": {
      ///       "option3": ...,
      ///       "option4": ...
      ///     }
      ///   }
      /// }
      /// \endcode
      ///
      ///
  */

  /*
    namespace {

      /// Checks:
      /// - no minimum orbit_branch_specs
      /// - warn for non-integer orbit_branch_specs
      /// - error if missing any in range [2, max(branch)]
      /// - for branch=1: ignore with warning
      /// - for branch 2+: max_length required,
      ///      max_length must be <= max_length for branch-1
      ///
      /// \code
      /// {
      ///   "cluster_specs": {
      ///     "method": "AperiodicClustersByMaxLength",
      ///     "kwargs": {
      ///       "standard": {
      ///         "orbit_branch_specs": {
      ///           "1": {
      ///             "max_length": 12.0,
      ///             "cutoff_radius": 6.0
      ///           },
      ///           "2": {
      ///             "max_length": 12.0,
      ///             "cutoff_radius": 6.0
      ///           },
      ///           "3": {
      ///             "max_length": 12.0,
      ///             "cutoff_radius": 6.0
      ///           }
      ///         }
      ///       },
      ///       "custom": [
      ///         {
      ///           "phenomenal": {
      ///             "coordinate_mode": "<>",
      ///             "sites": [...],
      ///             "equivalence": "prim", // "scel", "config"
      ///             "scelname": "<scelname>",
      ///             "configname": "<configname>"
      ///           },
      ///           "orbit_branch_specs": {}
      ///           "orbit_specs": [...(coordinates relative to phenom)...]
      ///         },
      ///         ...
      ///       ]
      ///     }
      ///   }
      /// }
      /// \endcode

      /// For parsing a phenomenal cluster
      ///
      /// \code
      /// {
      ///   "phenomenal": {
      ///     "coordinate_mode": "<>",
      ///     "sites": [...],
      ///     "equivalence": "prim", // one of ["prim", "scel", "config"]
      ///     "scelname": "<scelname>",  // required if "equivalence" == "scel"
      ///     "configname": "<configname>", // required if "equivalence" == "config"
      ///   }
      /// }
      /// \endcode
      template<typename PhenomenalType>
      struct AperiodicPhenomenalParser : KwargsParser {

        const PrimClex& primclex;
        PhenomenalType phenom;
        std::string equivalence;
        std::string scelname;
        std::string configname;

        AperiodicPhenomenalParser(
            const PrimClex& _primclex,
            jsonParser& _input,
            fs::path _path,
            bool _required):
          KwargsParser(_input, _path, _required),
          primclex(_primclex),
          phenom(_primclex.prim()) {

          // phenom:
          // - check equivalency wrt "prim" symmetry (prim.factor_group(), PrimSymCompare)
          // - check equivalency wrt "scel" symmetry (supercell.factor_group(), ScelSymCompare)
          // - check equivalency wrt "config" symmetry (config.factor_group(), ScelSymCompare)

          // generate custom local clusters:
          // OrbitGenerators<OrbitType> custom_generators(phenom, generating_sym_g, sym_compare)

          try {
            from_json(phenom, *self_it, primclex.crystallography_tol());
          }
          catch(std::exception& e) {
            std::stringstream msg;
            msg << "Error: could not read phenomenal " << traits<PhenomenalType>::name << ": "
                << e.what() << std::endl;
            error.insert(msg.str());
            continue;
          }

          require_enum<EQUIVALENCE_TYPE>();
        }
      };

      /// One of these for 'standard' specs, and one for each 'custom' specs (by phenomemal cluster)
      struct AperiodicOrbitBranchSpecsParser : OrbitBranchSpecsParser {

        AperiodicOrbitBranchSpecsParser(jsonParser& _input, fs::path _path, bool _required) :
          OrbitBranchSpecsParser(_input, _path, _required) {

          if(self_it != parent().end()) {
            for(auto it=self_it->begin(); it!=self_it->end(); ++it) {
              if(!warn_non_integer_branch(it)) {
                continue;
              }
              int branch = branch_to_int(it);
              if(branch < 1) { warn_unnecessary_branch(it); continue; }
              if(branch > 0) {
                require_option<double>(it, "max_length");
                require_option<double>(it, "cutoff_radius");
              }
              if(branch > 1) {
                require_nonincreasing<double>(it, "max_length");
                require_nonincreasing<double>(it, "cutoff_radius");
              }
            }
          }
        }
      };

      template<typename PhenomenalType, typename OrbitType>
      struct AperiodicOrbitSpecsParser : KwargsParser {

        const PrimClex& primclex;

        std::map<PhenomenalType, std::vector<typename OrbitType::Element> > phenom;
        typedef typename OrbitType::SymCompareType SymCompareType;

        AperiodicOrbitSpecsParser(
            const PrimClex& _primclex,
            const SymGroup& _generating_grp,
            const SymCompareType& _sym_compare,
            jsonParser& _input,
            fs::path _path,
            bool _required):
          KwargsParser(_input, _path, _required),
          primclex(_primclex),
          custom_generators(_generating_grp, _sym_compare) {

          // --- aperiodic clusters ---
          // read phenomenal cluster from specs
          // read custom generators from specs, referenced to phenomenal cluster location
          // check if subclusters should be included (yes by default)
          // generate custom generators

          // phenom:
          // - check equivalency wrt "prim" symmetry (prim.factor_group(), PrimSymCompare)
          // - check equivalency wrt "scel" symmetry (supercell.factor_group(), ScelSymCompare)
          // - check equivalency wrt "config" symmetry (config.factor_group(), ScelSymCompare)

          // generate custom local clusters:
          // OrbitGenerators<OrbitType> custom_generators(phenom, generating_sym_g, sym_compare)

          auto orbit_specs_it = parent().find(name());
          if(orbit_specs_it != parent().end()) {

            // for each custom orbit
            for(Index i=0; i<orbit_specs_it->size(); ++i) {
              const jsonParser& val = (*orbit_specs_it)[i];

              // read orbit generating cluster from specs
              typename OrbitType::Element input_cluster(primclex.prim());
              try {
                from_json(input_cluster, val, primclex.crystallography_tol());
              }
              catch(std::exception& e) {
                std::stringstream msg;
                msg << "Error in '" << name() << "'[" << i << "] reading cluster: "
                    << e.what() << std::endl << val << std::endl;
                error.insert(msg.str());
                continue;
              }

              try {
                // check if subclusters should be included (yes by default)
                auto f_it = val.find("include_subclusters");
                if(f_it == val.end() ||
                   (f_it != val.end() && f_it->get<bool>())) {
                  insert_subcluster_generators(input_cluster, custom_generators, primclex.log());
                }
                else {
                  custom_generators.insert(input_cluster);
                }
              }
              catch(std::exception& e) {
                std::stringstream msg;
                msg << "Error in '" << name() << "'[" << i << "] constructing generators: "
                    << e.what() << std::endl << val << std::endl;
                error.insert(msg.str());
                continue;
              }
            }
          }
        }
      };


      template<typename PhenomenalType, typename OrbitType>
      struct AperiodicClusterSpecsParser {

        //AperiodicClusterSpecsParserImpl<OrbitType> standard;
        //std::map<PhenomenalType, AperiodicClusterSpecsParserImpl<OrbitType>> custom;

      };

      template<typename PhenomenalType, typename OrbitType>
      using LocalClusterSpecsParser = AperiodicClusterSpecsParser<PhenomenalType, OrbitType>;
    }
    */

  /* -- Cluster Orbit generating function definitions ------------------------------------- */

  /// \brief Output the neighborhood of UnitCellCoord within max_radius of any site in unit cell
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


  /// \brief Output the neighborhood of DiffusionTransformation within max_radius of any site in transformation
  ///
  /// \param diff_trans DiffusionTransformation
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
    const Kinetics::DiffusionTransformation &diff_trans,
    double max_radius,
    std::function<bool (CoordType)> site_filter,
    OutputIterator result,
    double xtal_tol) {

    int max_low_shift = 0;
    int max_high_shift = 0;
    for(auto it = diff_trans.specie_traj().begin(); it != diff_trans.specie_traj().end(); it++) {
      Eigen::Vector3l vec = it->from.uccoord.unitcell();
      if(vec.maxCoeff() > max_high_shift) {
        max_high_shift = vec.maxCoeff();
      }
      if(vec.minCoeff() < max_low_shift) {
        max_low_shift = vec.minCoeff();
      }
    }

    auto dim = diff_trans.prim().lattice().enclose_sphere(max_radius);

    Eigen::Vector3i ones(1, 1, 1);
    EigenCounter<Eigen::Vector3i> grid_count(
      -dim + (max_low_shift * ones).cast<int>(),
      dim + (max_high_shift * ones).cast<int>(),
      Eigen::Vector3i::Constant(1));

    ///lattice scaling
    Coordinate lat_point(diff_trans.prim().lattice());

    const auto &basis = diff_trans.prim().basis;

    do {
      lat_point.frac() = grid_count().cast<double>();

      for(auto it = basis.begin(); it != basis.end(); ++it) {
        if(!site_filter(*it)) {
          continue;
        }

        Coordinate test(*it + lat_point);
        UnitCellCoord tmp(diff_trans.prim(), test, xtal_tol);
        if(dist_to_path(diff_trans, tmp) < max_radius && dist_to_path(diff_trans, tmp) > xtal_tol) {
          *result++ = UnitCellCoord(diff_trans.prim(), test, xtal_tol);
        }
        if(dist_to_path(diff_trans, tmp) <= xtal_tol) {
          auto spec_it = diff_trans.specie_traj().begin();
          for(; spec_it != diff_trans.specie_traj().end(); ++spec_it) {
            if(spec_it->to.uccoord == tmp || spec_it->from.uccoord == tmp) {
              break;
            }
          }
          if(spec_it == diff_trans.specie_traj().end()) {
            *result++ = UnitCellCoord(diff_trans.prim(), test, xtal_tol);
          }
        }
      }
    }
    while(++grid_count);
    return result;
  }


  /// \brief Check if the periodic image of a neighborhood overlaps itself
  ///
  /// \param local_orbits Vector of IntegralCluster orbits defining the local neighborhood
  /// \param scel A supercell
  ///
  template<typename OrbitType>
  bool has_local_neighborhood_overlap(std::vector<OrbitType> &local_orbits, const Supercell &scel) {
    std::set<int> present;
    std::set<UnitCellCoord> coords;
    for(auto &orbit : local_orbits) {
      for(auto &cluster : orbit) {
        for(int i = 0; i < cluster.size(); ++i) {
          coords.insert(cluster[i]);
        }
      }
    }
    for(auto &coord : coords) {
      if(!present.insert(scel.linear_index(coord)).second) {
        return true;
      }
    }
    //If no set insertion collision then no problems
    return false;
  }


  /// \brief Filter supercells to include only those with no local neighborhood overlap
  ///
  /// \param local_orbits Vector of IntegralCluster orbits defining the local neighborhood
  /// \param scel_options A vector of supercells
  ///
  template<typename OrbitType>
  std::vector<Supercell> viable_supercells(std::vector<OrbitType> &local_orbits, std::vector<Supercell> scel_options) {
    std::vector<Supercell> results;
    for(auto &scel : scel_options) {
      if(!has_local_neighborhood_overlap(local_orbits, scel)) {
        results.push_back(scel);
      }
    }
    return results;
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
    template<typename OrbitType>
    OrbitGenerators<OrbitType> &_insert_null_cluster_generator(
      const Structure &prim,
      OrbitGenerators<OrbitType> &generators) {

      typedef typename OrbitType::Element cluster_type;

      // create a prototype cluster
      cluster_type test(prim);
      generators.insert(test);

      return generators;
    }

    /// \brief Insert the generating elements for the asymmetric unit, including all sites
    ///
    /// \param prim A PrimType
    /// \param generators OrbitGenerators<OrbitType> instance collecting orbit
    ///        generating elements
    ///
    /// \relates IntegralCluster
    ///
    template<typename OrbitType>
    OrbitGenerators<OrbitType> &_insert_asymmetric_unit_generators(
      const Structure &prim,
      OrbitGenerators<OrbitType> &generators) {

      typedef typename OrbitType::Element cluster_type;

      // for all sites in the basis
      for(int i = 0; i < prim.basis.size(); i++) {
        // create a prototype cluster
        cluster_type test(prim);
        test.elements().push_back(UnitCellCoord(prim, i, UnitCell(0, 0, 0)));
        generators.insert(test);
      }

      return generators;
    }

    /// \brief Use orbits of size n to insert generating elements for orbits of size n+1
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
    template<typename OrbitType, typename OrbitGeneratorIterator>
    OrbitGenerators<OrbitType> &_insert_next_orbitbranch_generators(
      OrbitGeneratorIterator begin,
      OrbitGeneratorIterator end,
      const OrbitBranchSpecs<OrbitType> &specs,
      OrbitGenerators<OrbitType> &generators,
      std::ostream &status) {

      typedef typename OrbitType::Element cluster_type;
      typedef OrbitType orbit_type;

      const auto &sym_compare = specs.sym_compare();
      const auto &filter = specs.filter();
      const auto &g = specs.generating_group();

      // print status messages
      std::string clean(100, ' ');

      // contains a pair of iterators over candidate UnitCellCoord
      auto candidate_sites = specs.candidate_sites();
      Index orig_size = generators.elements.size();

      // for each orbit generator of size n
      for(auto orbit_generator_it = begin; orbit_generator_it != end; ++orbit_generator_it) {

        // print status messages
        status << clean << '\r'
               << "  Calculating orbit branch " << orbit_generator_it->size() + 1
               << ":  Expanding orbit " << std::distance(begin, orbit_generator_it)
               << " / " << std::distance(begin, end)
               << "  of branch " << orbit_generator_it->size()
               << ".  New orbits: " << generators.elements.size() - orig_size << std::flush;

        // by looping over each site in the grid,
        for(auto site_it = candidate_sites.first; site_it != candidate_sites.second; ++site_it) {

          // don't duplicate sites in cluster
          if(contains(*orbit_generator_it, *site_it)) {
            continue;
          }

          // create a test cluster from prototype
          cluster_type test(*orbit_generator_it);

          // add the new site
          test.elements().push_back(*site_it);
          // filter clusters
          if(!filter(test)) {
            continue;
          }

          // try inserting test (only uniques will be kept)
          generators.insert(test);
        }
      }

      // print status messages
      status << clean << '\r';

      return generators;
    }
  }


  // ------- Generating asymmetric unit orbits ---------------------------------

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
  OrbitOutputIterator make_next_orbitbranch(OrbitInputIterator begin,
                                            OrbitInputIterator end,
                                            const OrbitBranchSpecs<OrbitType> &specs,
                                            OrbitOutputIterator result,
                                            std::ostream &status) {

    /// Construct an OrbitGenerators object to collect orbit generating elements
    OrbitGenerators<OrbitType> generators(specs.generating_group(), specs.sym_compare());

    /// Use OrbitBranchSpecs to insert orbit generating elements for the next orbitbranch
    _insert_next_orbitbranch_generators(
      prototype_iterator(begin),
      prototype_iterator(end),
      specs,
      generators,
      status);

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
  template<typename OrbitBranchSpecsIterator, typename OrbitOutputIterator>
  OrbitOutputIterator make_orbits(
    OrbitBranchSpecsIterator begin,
    OrbitBranchSpecsIterator end,
    const std::vector<IntegralCluster> &custom_generators,
    OrbitOutputIterator result,
    std::ostream &status) {

    typedef typename OrbitOutputIterator::container_type container_type;
    typedef typename container_type::value_type orbit_type;

    // -- Collect orbit generating elements for each branch, to generate the next branch
    std::vector<OrbitGenerators<orbit_type> > generators;
    generators.reserve(std::distance(begin, end));

    // -- Combines all orbit branches
    OrbitGenerators<orbit_type> all_generators(begin->generating_group(), begin->sym_compare());

    // branches should be ordered already, so insert with end as hint
    auto insert_branch = [&](OrbitGenerators<orbit_type> &all, OrbitGenerators<orbit_type> &branch) {
      for(auto it = branch.elements.begin(); it != branch.elements.end(); ++it) {
        all.elements.insert(all.elements.end(), *it);
      }
    };
    // -- construct null cluster orbit

    auto specs_it = begin;
    if(specs_it != end) {

      generators.emplace_back(begin->generating_group(), begin->sym_compare());
      _insert_null_cluster_generator(specs_it->prim(), generators.back());
      ++specs_it;
      insert_branch(all_generators, generators.back());
    }
    // -- construct additional branches
    // print status messages
    std::string clean(100, ' ');

    // generate orbit branches 1+ using the previously generated branch:
    auto prev_gen = generators.begin();
    while(specs_it != end) {

      // print status message
      status << clean << '\r' << "Calculating orbit branch "
             << std::distance(begin, specs_it) << "\r" << std::flush;


      generators.emplace_back(specs_it->generating_group(), specs_it->sym_compare());
      _insert_next_orbitbranch_generators(
        prev_gen->elements.begin(),
        prev_gen->elements.end(),
        *specs_it,
        generators.back(),
        status);
      ++specs_it;
      ++prev_gen;
      insert_branch(all_generators, generators.back());
    }

    // -- add custom orbit generators

    for(int i = 0; i < custom_generators.size(); ++i) {
      status << clean << '\r' << "Adding custom orbit "
             << i << "/"
             << custom_generators.size() << "\r" << std::flush;

      all_generators.insert(custom_generators[i]);
    }

    // make orbits
    return all_generators.make_orbits(result);
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
  template<typename OrbitType>
  OrbitGenerators<OrbitType> &insert_subcluster_generators(
    typename OrbitType::Element cluster,
    OrbitGenerators<OrbitType> &generators,
    std::ostream &status) {

    typedef typename OrbitType::Element cluster_type;
    typedef OrbitType orbit_type;

    // Construct a functor that returns takes a cluster and returns it in a canonical form
    CanonicalGenerator<orbit_type> canonical_generator(generators.group, generators.sym_compare);

    // SubClusterGenerator allows iterating over subclusters (includes null and original cluster)
    SubClusterGenerator<cluster_type> it(cluster);
    SubClusterGenerator<cluster_type> end;

    while(it != end) {
      generators.elements.insert(canonical_generator(*it++));
    }

    return generators;
  }

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
    const Structure &prim,
    const std::function<bool (Site)> &site_filter,
    double xtal_tol,
    OrbitOutputIterator result,
    std::ostream &status) {

    typedef PrimPeriodicOrbit<IntegralCluster> orbit_type;
    typedef typename orbit_type::Element cluster_type;

    const SymGroup &generating_grp = prim.factor_group();
    PrimPeriodicSymCompare<IntegralCluster> sym_compare(xtal_tol);

    std::vector<UnitCellCoord> candidate_sites;
    for(int i = 0; i < prim.basis.size(); ++i) {
      if(site_filter(prim.basis[i])) {
        candidate_sites.emplace_back(prim, i, 0, 0, 0);
      }
    }

    OrbitBranchSpecs<orbit_type> specs(
      prim,
      candidate_sites.begin(),
      candidate_sites.end(),
      generating_grp,
    [](const cluster_type & test) {
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
    const Structure &prim,
    const std::vector<double> &max_length,
    const std::vector<IntegralCluster> &custom_generators,
    const std::function<bool (Site)> &site_filter,
    double xtal_tol,
    OrbitOutputIterator result,
    std::ostream &status) {

    typedef PrimPeriodicOrbit<IntegralCluster> orbit_type;
    typedef typename orbit_type::Element cluster_type;

    const SymGroup &generating_grp = prim.factor_group();
    PrimPeriodicSymCompare<IntegralCluster> sym_compare(xtal_tol);

    // collect OrbitBranchSpecs here
    std::vector<OrbitBranchSpecs<orbit_type> > specs;

    // collect the environment of sites here
    std::vector<UnitCellCoord> candidate_sites;

    // --- add specs for null cluster orbit ------------------
    if(max_length.size() >= 1) {
      specs.emplace_back(prim,
                         candidate_sites.begin(),
                         candidate_sites.end(),
                         generating_grp,
      [](const cluster_type & test) {
        return true;
      },
      sym_compare);
    }


    // --- add specs for asymmetric unit orbit ------------------
    if(max_length.size() >= 2) {
      for(int i = 0; i < prim.basis.size(); ++i) {
        if(site_filter(prim.basis[i])) {
          candidate_sites.emplace_back(prim, i, 0, 0, 0);
        }
      }
      specs.emplace_back(prim,
                         candidate_sites.begin(),
                         candidate_sites.end(),
                         generating_grp,
      [](const cluster_type & test) {
        return true;
      },
      sym_compare);
    }

    // --- add specs for additional orbit branches ------------------
    for(auto it = max_length.begin() + 2; it != max_length.end(); ++it) {

      // construct the neighborhood of sites to consider for the orbit
      candidate_sites.clear();
      neighborhood(prim, *it, site_filter, std::back_inserter(candidate_sites), xtal_tol);

      auto max_length_filter = [ = ](const cluster_type & test) {
        return test.invariants().displacement().back() < *it;
      };

      specs.emplace_back(prim,
                         candidate_sites.begin(),
                         candidate_sites.end(),
                         generating_grp,
                         max_length_filter,
                         sym_compare);
    }

    // now generate orbits
    return make_orbits(specs.begin(), specs.end(), custom_generators, result, status);

  }

  /// \brief Generate Orbit<IntegralCluster> around DiffusionTransformation
  /// from bspecs.json-type JSON input file
  ///
  /// \param prim Primitive structure
  /// \param bspecs jsonParser containing bspecs.json contents
  /// \param site_filter A filter function that returns true for Site that
  ///        should be considered for the neighborhood (i.e. to check the number
  ///        of components)
  /// \param xtal_tol Crystallography tolerance
  /// \param result An output iterator for Orbit
  /// \param status Stream for status messages
  ///
  /// - Uses prim.factor_group as the generating group
  /// - Uses PrimPerioidicSymCompare<IntegralCluster>(xtal_tol) for cluster equivalence
  /// - Converts input to max_length and custom_generators and calls make_orbits
  ///
  /// \relates IntegralCluster
  template<typename OrbitOutputIterator>
  OrbitOutputIterator make_prim_periodic_orbits(
    const Structure &prim,
    const jsonParser &bspecs,
    const std::function<bool (Site)> &site_filter,
    double xtal_tol,
    OrbitOutputIterator result,
    std::ostream &status) {

    typedef PrimPeriodicOrbit<IntegralCluster> orbit_type;
    typedef typename orbit_type::Element cluster_type;

    // read max_length from bspecs
    std::vector<double> max_length = max_length_from_bspecs(bspecs);

    // collect custom orbit generating clusters in 'generators'
    PrimPeriodicSymCompare<IntegralCluster> sym_compare(xtal_tol);
    OrbitGenerators<orbit_type> generators(prim.factor_group(), sym_compare);

    if(bspecs.contains("orbit_specs")) {

      // for each custom orbit
      for(auto it = bspecs["orbit_specs"].begin(); it != bspecs["orbit_specs"].end(); ++it) {

        // read orbit generating cluster from bspecs
        cluster_type input_cluster(prim);
        from_json(input_cluster, *it, xtal_tol);

        // check if subclusters should be included (yes by default)
        auto f_it = it->find("include_subclusters");
        if(f_it == it->end() ||
           (f_it != it->end() && f_it->get<bool>())) {
          insert_subcluster_generators(input_cluster, generators, status);
        }
        else {
          generators.insert(input_cluster);
        }
      }
    }

    std::vector<cluster_type> custom_generators(generators.elements.begin(), generators.elements.end());

    return make_prim_periodic_orbits(prim, max_length, custom_generators, site_filter, xtal_tol, result, status);

  }

  /// \brief Generate Orbit<IntegralCluster> by specifying max cluster length for each branch
  /// by specifying max cluster length for each branch and cut off radius for local environment
  ///
  /// \param diff_trans DiffusionTransformation
  /// \param generating_group Invariant group of diff_trans used to generate local orbits
  /// \param sym_compare Comparison used to check cluster equivalence
  /// \param cutoff_radius max radius from the transformation that defines the local environment
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
  /// - Figures out candidate_sites from max_length and site_filter input to
  ///   create OrbitBranchSpecs and calls make_orbits
  ///
  /// \relates IntegralCluster
  template<typename OrbitOutputIterator, typename SymCompareType>
  OrbitOutputIterator make_local_orbits(
    const Kinetics::DiffusionTransformation &diff_trans,
    const SymGroup &generating_group,
    const SymCompareType &sym_compare,
    const std::vector<double> &cutoff_radius,
    const std::vector<double> &max_length,
    const std::vector<IntegralCluster> &custom_generators,
    const std::function<bool (Site)> &site_filter,
    double xtal_tol,
    OrbitOutputIterator result,
    std::ostream &status) {

    typedef typename OrbitOutputIterator::container_type container_type;
    typedef typename container_type::value_type orbit_type;
    typedef typename orbit_type::Element cluster_type;

    // collect OrbitBranchSpecs here
    std::vector<OrbitBranchSpecs<orbit_type> > specs;

    // collect the environment of sites here
    std::vector<UnitCellCoord> candidate_sites;

    // --- add specs for null cluster orbit ------------------
    if(max_length.size() >= 1) {
      specs.emplace_back(diff_trans.prim(),
                         candidate_sites.begin(),
                         candidate_sites.end(),
                         generating_group,
      [](const cluster_type & test) {
        return true;
      },
      sym_compare);
    }
    // --- add specs for asymmetric unit orbit ------------------
    if(max_length.size() >= 2) {
      neighborhood(diff_trans, cutoff_radius[1], site_filter, std::back_inserter(candidate_sites), xtal_tol);

      specs.emplace_back(diff_trans.prim(),
                         candidate_sites.begin(),
                         candidate_sites.end(),
                         generating_group,
      [](const cluster_type & test) {
        return true;
      },
      sym_compare);
    }
    int idx = 1;
    // --- add specs for additional orbit branches ------------------
    for(auto it = max_length.begin() + 2; it != max_length.end(); ++it) {
      ++idx;
      // construct the neighborhood of sites to consider for the orbit
      candidate_sites.clear();
      neighborhood(diff_trans, cutoff_radius[idx], site_filter, std::back_inserter(candidate_sites), xtal_tol);

      auto max_length_filter = [ = ](const cluster_type & test) {
        double check = test.invariants().displacement().back() ;
        return check < *it;
      };

      specs.emplace_back(diff_trans.prim(),
                         candidate_sites.begin(),
                         candidate_sites.end(),
                         generating_group,
                         max_length_filter,
                         sym_compare);
    }

    // now generate orbits
    return make_orbits(specs.begin(), specs.end(), custom_generators, result, status);

  }


  /// \brief Generate Orbit<IntegralCluster> from bspecs.json-type JSON input file
  ///
  /// \param diff_trans Kinetics::DiffusionTransformation
  /// \param generating_group Invariant group of diff_trans used to generate local orbits
  /// \param sym_compare Comparison used to check cluster equivalence
  /// \param bspecs jsonParser containing bspecs.json contents
  /// \param site_filter A filter function that returns true for Site that
  ///        should be considered for the neighborhood (i.e. to check the number
  ///        of components)
  /// \param xtal_tol Crystallography tolerance
  /// \param result An output iterator for Orbit
  /// \param status Stream for status messages
  ///
  /// - Reads bspecs input for cutoff_radius, max_length, and custom_generators and calls
  ///   make_local_orbits overload.
  ///
  /// \relates IntegralCluster
  template<typename OrbitOutputIterator, typename SymCompareType>
  OrbitOutputIterator make_local_orbits(
    const Kinetics::DiffusionTransformation &diff_trans,
    const SymGroup &generating_group,
    const SymCompareType &sym_compare,
    const jsonParser &bspecs,
    const std::function<bool (Site)> &site_filter,
    double xtal_tol,
    OrbitOutputIterator result,
    std::ostream &status) {

    typedef typename OrbitOutputIterator::container_type container_type;
    typedef typename container_type::value_type orbit_type;
    typedef typename orbit_type::Element cluster_type;

    std::vector<double> max_length;
    std::vector<double> cutoff_radius;
    try {
      // read max_length from bspecs
      max_length = max_length_from_bspecs(bspecs);
      // read cutoff_radius from bspecs
      cutoff_radius = cutoff_radius_from_bspecs(bspecs);
    }
    catch(std::exception &e) {
      default_err_log().error("In make_local_orbits (from bspecs)");
      default_err_log() << "Error reading max_length or cutoff_radius" << std::endl;
      default_err_log() << e.what() << std::endl;
      throw e;
    }

    // collect custom orbit generating clusters in 'generators'
    OrbitGenerators<orbit_type> generators(generating_group, sym_compare);

    try {
      if(bspecs.contains("orbit_specs")) {

        // for each custom orbit
        for(auto it = bspecs["orbit_specs"].begin(); it != bspecs["orbit_specs"].end(); ++it) {

          // read orbit generating cluster from bspecs
          cluster_type input_cluster(diff_trans.prim());
          from_json(input_cluster, *it, xtal_tol);

          // check if subclusters should be included (yes by default)
          auto f_it = it->find("include_subclusters");
          if(f_it == it->end() ||
             (f_it != it->end() && f_it->get<bool>())) {
            insert_subcluster_generators(input_cluster, generators, status);
          }
          else {
            generators.insert(input_cluster);
          }
        }
      }
    }
    catch(std::exception &e) {
      default_err_log().error("In make_local_orbits (from bspecs)");
      default_err_log() << "Error getting custom orbit generators" << std::endl;
      default_err_log() << e.what() << std::endl;
      throw e;
    }

    std::vector<cluster_type> custom_generators(generators.elements.begin(), generators.elements.end());

    return make_local_orbits(diff_trans, generating_group, sym_compare, cutoff_radius, max_length, custom_generators, site_filter, xtal_tol, result, status);

  }


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
