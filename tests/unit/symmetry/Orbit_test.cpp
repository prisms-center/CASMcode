#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/symmetry/InvariantSubgroup.hh"
#include "casm/symmetry/InvariantSubgroup_impl.hh"
#include "casm/symmetry/Orbit.hh"
#include "casm/symmetry/Orbit_impl.hh"

/// What is being used to test it:
#include "Common.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"

using namespace CASM;

namespace test {
  /// \brief Returns vector containing sorted indices of SymOp in first column of equivalence_map
  ///
  /// - This is what the first column of equivalence_map will look like if Element
  /// 'proto' is the prototype.
  /// - Uses index into generating_group
  template<typename Element, typename SymCompareType, typename EquivContainer>
  std::vector<Index> _sorter(
    const Element &proto,
    const EquivContainer &equiv,
    const SymGroup &g,
    const SymCompareType &sym_compare) {

    auto equal = [&](const Element & A, const Element & B) {
      return sym_compare.equal(A, B);
    };

    int count = 0;
    std::vector<Index> sorter(equiv.size(), -1);

    for(Index op_index = 0; op_index != g.size(); ++op_index) {
      Index i = find_index(equiv, sym_compare.prepare(copy_apply(g[op_index], proto)), equal);
      if(sorter[i] == -1) {
        sorter[i] = op_index;
        count++;
        if(count == equiv.size()) {
          std::sort(sorter.begin(), sorter.end());
          return sorter;
        }
      }
    }
    throw std::runtime_error("Error generating equivalence map");
  }
}

BOOST_AUTO_TEST_SUITE(OrbitTest)

BOOST_AUTO_TEST_CASE(Test0) {
  test::ZrOProj proj;
  proj.check_init();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);
  const Structure &prim = primclex.prim();
  const Lattice &lat = prim.lattice();
  Supercell prim_scel(&primclex, Eigen::Matrix3i::Identity());

  BOOST_CHECK_EQUAL(true, true);

  // Make PrimPeriodicIntegralClusterOrbit
  fs::path bspecs_path = "tests/unit/kinetics/bspecs_0.json";
  jsonParser bspecs {bspecs_path};

  std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
  make_prim_periodic_orbits(
    primclex.prim(),
    bspecs,
    alloy_sites_filter,
    primclex.crystallography_tol(),
    std::back_inserter(orbits),
    primclex.log());
  BOOST_CHECK_EQUAL(true, true);

  // Make cluster groups & check size, based on reduced symmetry of a vol 2 Supercell
  {
    // Make vol 2 supercell & and background configuration
    Eigen::Vector3d a, b, c;
    std::tie(a, b, c) = primclex.prim().lattice().vectors();
    Supercell scel_vol2 {&primclex, Lattice(2 * a, 1 * b, 1 * c)};
    Configuration config(scel_vol2);

    // Scel sym_compare
    ScelPeriodicSymCompare<IntegralCluster> scel_sym_compare(
      scel_vol2.prim_grid(),
      scel_vol2.crystallography_tol());


    // Get the config factor group (should just be the Supercell factor group)
    std::vector<PermuteIterator> _config_fg = config.factor_group();
    SymGroup config_fg = make_sym_group(_config_fg.begin(), _config_fg.end());

    // - these values have not been checked for correctness, they just check for consistency
    std::vector<Index> expected_cluster_group_size = {16, 4, 4, 4, 2, 4, 2, 2, 1,
                                                      4, 2, 2, 2, 4, 2, 2, 2, 1, 2, 1, 1, 2, 4, 2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1,
                                                      1, 1, 2, 1, 1, 1, 1, 4, 2, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1,
                                                      2, 4, 2, 1, 2, 2, 2, 1, 2, 1, 1, 1, 2, 2, 1
                                                     };
    std::vector<Index> cluster_group_size;
    Index index = 0;
    for(const auto &orbit : orbits) {

      //std::cout << "\n ---------- Orbit: " << index << " ----------- \n" << std::endl;

      // Test make_invariant_subgroup using orbit generators
      {
        SymGroup cluster_group = make_invariant_subgroup(
                                   orbit.prototype(),
                                   config_fg,
                                   scel_sym_compare);
        cluster_group_size.push_back(cluster_group.size());
        BOOST_CHECK_EQUAL(expected_cluster_group_size[index], cluster_group.size());
      }

      // Test the components of the orbit constructor
      {
        // data
        const auto &generating_element = orbit.prototype();
        const auto &g = config_fg;
        const auto &m_sym_compare = scel_sym_compare;
        typedef IntegralCluster Element;
        std::vector<Element> m_element;
        multivector<SymOp>::X<2> m_equivalence_map;

        BOOST_CHECK_EQUAL(g.is_group(), true);

        // --- recreate steps in Orbit constructor

        // define functions
        BOOST_CHECK_EQUAL(true, true);
        auto prepare = [&](const Element & A) {
          return m_sym_compare.prepare(A);
        };
        auto compare = [&](const Element & A, const Element & B) {
          return m_sym_compare.compare(A, B);
        };
        auto equal = [&](const Element & A, const Element & B) {
          return m_sym_compare.equal(A, B);
        };
        BOOST_CHECK_EQUAL(true, true);

        // generate equivalents
        std::set<Element, decltype(compare)> t_equiv(compare);
        for(const auto &op : g) {
          t_equiv.insert(prepare(copy_apply(op, generating_element)));
        }
        BOOST_CHECK_EQUAL(true, true);

        // sort element using each element's first equivalence map column to find prototype
        std::set<std::pair<Element, std::vector<Index> >, Orbit_impl::_EqMapCompare<Element> > _set;
        for(const auto &e : t_equiv) {
          _set.insert(std::make_pair(e, Orbit_impl::_sorter(e, t_equiv, g, m_sym_compare)));
        }
        BOOST_CHECK_EQUAL(true, true);

        // use _set.begin()->first for prototype, use _set.begin()->second to generate equiv
        for(auto op_index : _set.begin()->second) {
          m_element.push_back(prepare(copy_apply(g[op_index], _set.begin()->first)));
        }
        BOOST_CHECK_EQUAL(true, true);
        BOOST_CHECK_EQUAL(m_element.size(), t_equiv.size());

        // generate equivalence map
        m_equivalence_map.resize(m_element.size());
        for(const auto &op : g) {
          Index i = find_index(m_element, prepare(copy_apply(op, m_element[0])), equal);
          m_equivalence_map[i].push_back(op);
        }
        BOOST_CHECK_EQUAL(true, true);


        // --- use actual Constructor ---
        BOOST_CHECK_EQUAL(true, true);
        Orbit<IntegralCluster, ScelPeriodicSymCompare<IntegralCluster>> suborbit(
                                                                       generating_element, g, m_sym_compare);
        BOOST_CHECK_EQUAL(true, true);
      }

      index++;
    }

    //test::print_computed_result(std::cout, "cluster_group_size", cluster_group_size);
    BOOST_CHECK_EQUAL(true, true);
  }
}

BOOST_AUTO_TEST_SUITE_END()
