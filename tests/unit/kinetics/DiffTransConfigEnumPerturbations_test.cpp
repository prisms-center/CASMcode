#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/kinetics/DiffTransConfigEnumPerturbations.hh"
#include "casm/kinetics/DiffTransEnumEquivalents.hh"

/// What is being used to test it:
#include "casm/clex/PrimClex.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/AppIO_impl.hh"
#include "Common.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/kinetics/DoFTransformation.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffusionTransformationEnum.hh"
#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/kinetics/SubOrbitGenerators.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/Orbit.hh"
#include "casm/casm_io/VaspIO.hh"

using namespace CASM;
using namespace test;

typedef Orbit <
Kinetics::DiffusionTransformation,
         Kinetics::PrimPeriodicDiffTransSymCompare > PrimPeriodicDiffTransOrbit;


//generating group: G {g0, g1, ...}
//orbit: O {e0, e1, ...} = G * proto
//sub-group of G: H {h0, h1, ...}

// ex. config.primitive() is cubic, has G && config is tetragonal, has H

// Orbit_g(proto, G) ->
//   Orbit_h0(C00'*proto, H), Orbit_h1(C10'*proto, H), ...

// all the unique cosets of H:
//   C0, C1, C2, ... have no elements in common, have all elements in G

// Ci0*(Sj*proto) may be equal to Cj0*(Sk*proto), where S is in invariant subgroup of proto

// so we want elements of G, g, such that g >= h*g*s, for all h in H, s in S



template<typename OrbitType>
SymGroup make_invariant_subgroup(const OrbitType &orbit, Index element_index = 0) {
  SymGroup result;
  result.set_lattice(orbit.prototype().lattice());
  const auto &map = orbit.equivalence_map();
  for(Index i = 0; i < orbit.equivalence_map()[0].size(); ++i) {
    result.push_back(map[0][i]*map[element_index][0]);
  }
  return result;
}

/// Implementation
template<typename Element, typename ElementOutputIterator>
ElementOutputIterator make_suborbit_generators(
  const Element &element,
  const SymGroup &invariant_subgroup,
  const SymGroup &group,
  const SymGroup &subgroup,
  ElementOutputIterator result) {

  // find "max" SymOp in each coset of the subgroup,
  //   excluding those cosets that cause duplicate generating elements
  for(Index i = 0; i < group.size(); ++i) {
    const SymOp &test_op = group[i];
    auto lambda = [&](const SymOp & op) {
      for(const auto &el_op : invariant_subgroup) {
        if(test_op.index() < (op * test_op * el_op).index()) {
          return true;
        }
      }
      return false;
    };
    // if test_op is max
    if(std::none_of(subgroup.begin(), subgroup.end(), lambda)) {
      // apply to element and construct suborbit generator
      *result++ = copy_apply(test_op, element);
    }
  }
  return result;
}

template<typename Element, typename SymCompareType, typename ElementOutputIterator>
ElementOutputIterator make_suborbit_generators(
  const Element &element,
  const SymCompareType &sym_compare,
  const SymGroup &group,
  const SymGroup &subgroup,
  ElementOutputIterator result) {

  SymGroup invariant_subgroup = make_invariant_subgroup(element, group, sym_compare);
  return make_suborbit_generators(element, invariant_subgroup, group, subgroup, result);
}

template<typename OrbitType, typename ElementOutputIterator>
ElementOutputIterator make_suborbit_generators(
  const OrbitType &orbit,
  const SymGroup &group,
  const SymGroup &subgroup,
  ElementOutputIterator result) {

  SymGroup invariant_subgroup = make_invariant_subgroup(orbit);
  return make_suborbit_generators(orbit.prototype(), invariant_subgroup, group, subgroup, result);
}

template<typename Element>
std::vector<PermuteIterator> make_invariant_subgroup(const Element &element, const Supercell &scel) {

  ScelPeriodicSymCompare<Element> sym_compare(scel.prim_grid(), scel.primclex().crystallography_tol());
  Element e(sym_compare.prepare(element));
  std::vector<PermuteIterator> result;
  auto it = scel.permute_begin();
  auto end = scel.permute_end();
  while(it != end) {
    auto test = sym_compare.prepare(copy_apply(it.sym_op(), e));
    if(sym_compare.equal(test, e)) {
      auto trans_it = scel.permute_it(0, scel.prim_grid().find(sym_compare.integral_tau()));
      result.push_back(trans_it * it);
    }
    ++it;
  }
  return result;
}

template<typename Element, typename ElementOutputIterator, typename PermuteIteratorIt>
ElementOutputIterator make_suborbit_generators(
  const Element &element,
  const Supercell &scel,
  PermuteIteratorIt subgroup_begin,
  PermuteIteratorIt subgroup_end,
  ElementOutputIterator result) {

  const PrimGrid &prim_grid = scel.prim_grid();

  std::vector<PermuteIterator> scel_inv_group = make_invariant_subgroup(element, scel);

  // find "max" permute in each coset of the subgroup [begin, end),
  //   excluding those cosets that cause duplicate generating elements
  auto test_it = scel.permute_begin();
  auto test_end = scel.permute_end();
  for(; test_it != test_end; ++test_it) {
    auto lambda = [&](const PermuteIterator & permute_it) {
      for(auto it = scel_inv_group.begin(); it != scel_inv_group.end(); ++it) {
        if(test_it < permute_it * (test_it * (*it))) {
          return true;
        }
      }
      return false;
    };

    // if test_it is max
    if(std::none_of(subgroup_begin, subgroup_end, lambda)) {
      // apply to prototype and construct suborbit
      *result++ = copy_apply(test_it.sym_op(), element);
    }
  }

  return result;
}

/// Prim-periodic generating elements -> Scel-periodic orbit generating elements in a Configuration
///
/// Uses config.supercell() PermuteIterators
template<typename OrbitIterator, typename ElementOutputIterator>
ElementOutputIterator make_suborbit_generators_slow(
  OrbitIterator begin,
  OrbitIterator end,
  const Configuration &config,
  ElementOutputIterator result) {

  typedef typename std::iterator_traits<OrbitIterator>::value_type OrbitType;
  typedef typename OrbitType::Element Element;

  const Structure &prim = config.prim();
  const Supercell &scel = config.supercell();
  std::vector<PermuteIterator> config_fg = config.factor_group();

  for(auto it = begin; it != end; ++it) {

    // get generating elements for prim->supercell orbit splitting
    std::vector<Element> scel_suborbit_generators;

    make_suborbit_generators(
      *it,
      prim.factor_group(),
      scel.factor_group(),
      std::back_inserter(scel_suborbit_generators));

    // get generating elements for supercell->config orbit splitting
    for(const auto &el : scel_suborbit_generators) {

      result = make_suborbit_generators(
                 el,
                 scel,
                 config_fg.begin(),
                 config_fg.end(),
                 result);
    }
  }
  return result;
}

/// Prim-periodic generating element -> Scel-periodic generating elements in a Configuration
///
/// Uses config.primitive().supercell() PermuteIterators, then splits again
/// for case that config has lower symmetry than config.primitive()
template<typename OrbitIterator, typename ElementOutputIterator>
ElementOutputIterator make_suborbit_generators(
  OrbitIterator begin,
  OrbitIterator end,
  const Configuration &config,
  ElementOutputIterator result) {

  typedef typename std::iterator_traits<OrbitIterator>::value_type OrbitType;
  typedef typename OrbitType::Element Element;

  Configuration prim_config = config.primitive().in_canonical_supercell();
  ScelPeriodicSymCompare<Element> sym_compare(
    config.supercell().prim_grid(),
    config.crystallography_tol());

  std::vector<Element> prim_config_suborbit_generators;

  make_suborbit_generators_slow(
    begin,
    end,
    prim_config,
    std::back_inserter(prim_config_suborbit_generators));

  for(const auto &el : prim_config_suborbit_generators) {
    result = make_suborbit_generators(
               el,
               sym_compare,
               prim_config.supercell().factor_group(),
               config.supercell().factor_group(),
               result);
  }

  return result;
}

// Test 1: PrimPeriodic -> ScelPeriodic orbits
template<typename OrbitType, typename ElementType>
void test_1(
  const Configuration &config,
  const Configuration &prim_config,
  const std::vector<OrbitType> &orbits,
  std::vector<ElementType> &scel_generators,
  std::vector<Index> &orbit_index,
  std::vector<Index> &scel_suborbit_size) {

  // useful SymGroups
  const Structure &prim = config.prim();
  const auto &prim_fg = prim.factor_group();
  const auto &config_scel_fg = config.supercell().factor_group();
  const auto &prim_config_scel_fg = prim_config.supercell().factor_group();

  std::cout << "Symmetry:\n";
  std::cout << "prim_fg.size(): " << prim_fg.size() << std::endl;
  std::cout << "config_scel_fg.size(): " << config_scel_fg.size() << std::endl;
  std::cout << "prim_config_scel_fg.size(): " << prim_config_scel_fg.size() << std::endl;

  std::cout << "\n!!! begin test 1" << std::endl;
  // Split IntegralCluster orbits according to config.supercell().factor_group()
  Index orbit_i = 0;
  for(const auto &orbit : orbits) {
    std::cout << "\n ----------------- \n";
    std::cout << "begin orbit " << orbit_i << "/" << orbits.size() << std::endl;
    std::vector<ElementType> suborbit_generators;
    make_suborbit_generators(
      orbit,
      prim_fg,
      prim_config_scel_fg,
      std::back_inserter(suborbit_generators));

    // Check results:
    std::cout << "Prototype, orbit " << orbit_i << "  prim-orbit.size(): "
              << orbit.size() << std::endl;
    //std::cout << orbit.prototype() << std::endl;
    std::cout << "Sub-orbit generators: " << std::endl;
    Index suborbit_i = 0;
    Index suborbit_size_sum = 0;
    for(const auto &el : suborbit_generators) {
      Orbit<ElementType, PrimPeriodicSymCompare<ElementType>> suborbit(el, prim_config_scel_fg, orbit.sym_compare());
      std::cout << "sub-orbit " << suborbit_i << "/" << suborbit_generators.size()
                << ", size: " << suborbit.size() << ":" << std::endl;
      //std::cout << el << std::endl;
      scel_suborbit_size.push_back(suborbit.size());
      orbit_index.push_back(orbit_i);
      suborbit_size_sum += suborbit.size();
      ++suborbit_i;
    }
    std::cout << "sum: " << suborbit_size_sum
              << "  prim-orbit.size(): "
              << orbit.size() << std::endl;

    BOOST_CHECK_EQUAL(orbit.size(), suborbit_size_sum);
    ++orbit_i;
    std::copy(suborbit_generators.begin(), suborbit_generators.end(), std::back_inserter(scel_generators));
  }

}

// Test 2: ScelPeriodic orbits -> Config orbits
template<typename OrbitType, typename ElementType>
void test_2(
  const Configuration &prim_config,
  const std::vector<OrbitType> &orbits,
  const std::vector<ElementType> &scel_generators,
  std::vector<ElementType> &config_generators,
  std::vector<Index> &orbit_index,
  std::vector<Index> &scel_suborbit_size) {

  const Structure &prim = prim_config.prim();
  std::vector<PermuteIterator> prim_config_fg = prim_config.factor_group();
  SymGroup _prim_config_fg = make_sym_group(prim.lattice(), prim_config_fg);
  ScelPeriodicSymCompare<ElementType> prim_config_scel_sym_compare(
    prim_config.supercell().prim_grid(),
    prim_config.crystallography_tol());

  std::vector<Index> config_suborbit_size;
  std::cout << "\n!!! begin test 2" << std::endl;
  // Split IntegralCluster orbits according to config occupation
  for(Index el_i = 0; el_i < scel_generators.size(); ++el_i) {
    std::cout << "\n ----------------- \n";
    std::cout << "begin scel_generator " << el_i << "/" << scel_generators.size() << std::endl;
    const auto &el = scel_generators[el_i];

    std::vector<ElementType> suborbit_generators;
    make_suborbit_generators(
      el,
      prim_config.supercell(),
      prim_config_fg.begin(),
      prim_config_fg.end(),
      std::back_inserter(suborbit_generators));

    // Check results:
    std::cout << "Generator, " << el_i << "  orbit, " << orbit_index[el_i] << "  scel-orbit.size() / prim-orbit.size(): "
              << scel_suborbit_size[el_i] << "/" << orbits[orbit_index[el_i]].size() << std::endl;
    //std::cout << el << std::endl;
    std::cout << "Sub-orbit generators: " << std::endl;
    Index config_suborbit_i = 0;
    Index config_suborbit_size_sum = 0;
    for(const auto &config_el : suborbit_generators) {
      Orbit<ElementType, ScelPeriodicSymCompare<ElementType>> suborbit(config_el, _prim_config_fg, prim_config_scel_sym_compare);
      std::cout << "config sub-orbit " << config_suborbit_i << "/" << suborbit_generators.size()
                << ", size: " << suborbit.size() << ":" << std::endl;
      //std::cout << config_el << std::endl;
      config_suborbit_size.push_back(suborbit.size());
      config_suborbit_size_sum += suborbit.size();
      ++config_suborbit_i;
    }
    std::cout << "sum: " << config_suborbit_size_sum
              << "  sub-orbit.size() * scel volume: "
              << scel_suborbit_size[el_i]*prim_config.supercell().volume() << std::endl;

    if(el.size() == 0) {
      BOOST_CHECK_EQUAL(1, config_suborbit_size_sum);
    }
    else {
      BOOST_CHECK_EQUAL(scel_suborbit_size[el_i]*prim_config.supercell().volume(), config_suborbit_size_sum);
    }
    std::copy(suborbit_generators.begin(), suborbit_generators.end(), std::back_inserter(config_generators));
  }
  std::cout << "  config_generators.size(): " << config_generators.size() << std::endl;
}

// Test 3: PrimPeriodic -> Config orbits,
//   using make_suborbit_generators_slow and primitive Configuration
template<typename OrbitType, typename ElementType>
void test_3(
  const Configuration &prim_config,
  const std::vector<OrbitType> &orbits,
  const std::vector<ElementType> &config_generators) {

  std::cout << "\n!!! begin test 3" << std::endl;
  std::vector<ElementType> config_generators_slow;
  {
    make_suborbit_generators_slow(
      orbits.begin(),
      orbits.end(),
      prim_config,
      std::back_inserter(config_generators_slow));
    std::cout << "  config_generators_slow.size(): " << config_generators_slow.size() << std::endl;
    BOOST_CHECK_EQUAL(config_generators.size(), config_generators_slow.size());
  }

}

// Test 4: PrimPeriodic -> Config orbits,
//   using make_suborbit_generators and non-primitive Configuration
//   (Need a new test config: current prim_config is lower symmetry than config)
template<typename OrbitType>
void test_4(
  const Configuration &config,
  const std::vector<OrbitType> &orbits) {

  typedef typename OrbitType::Element ElementType;

  std::cout << "\n!!! begin test 4" << std::endl;
  std::vector<ElementType> config_generators_nonprim_slow;
  std::vector<ElementType> config_generators_nonprim;
  {
    make_suborbit_generators_slow(
      orbits.begin(),
      orbits.end(),
      config,
      std::back_inserter(config_generators_nonprim_slow));
    std::cout << "  config_generators_nonprim_slow.size(): " << config_generators_nonprim_slow.size() << std::endl;

    make_suborbit_generators(
      orbits.begin(),
      orbits.end(),
      config,
      std::back_inserter(config_generators_nonprim));
    std::cout << "  config_generators_nonprim.size(): " << config_generators_nonprim.size() << std::endl;

    BOOST_CHECK_EQUAL(config_generators_nonprim_slow.size(), config_generators_nonprim.size());
  }
}



BOOST_AUTO_TEST_SUITE(DiffTransConfigEnumPerturbationsTest)

BOOST_AUTO_TEST_CASE(Test0) {

  test::ZrOProj proj;
  proj.check_init();
  proj.check_composition();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);
  const Structure &prim = primclex.prim();
  const Lattice &lat = prim.lattice();

  fs::path bspecs_path = "tests/unit/kinetics/bspecs_0.json";
  jsonParser bspecs {bspecs_path};

  // Make PrimPeriodicIntegralClusterOrbit
  std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
  make_prim_periodic_orbits(
    primclex.prim(),
    bspecs,
    alloy_sites_filter,
    primclex.crystallography_tol(),
    std::back_inserter(orbits),
    primclex.log());

  print_clust(orbits.begin(), orbits.end(), std::cout, PrototypePrinter<IntegralCluster>());

  // Make PrimPeriodicDiffTransOrbit
  std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
  Kinetics::make_prim_periodic_diff_trans_orbits(
    orbits.begin() + 2,
    orbits.begin() + 4,
    primclex.crystallography_tol(),
    std::back_inserter(diff_trans_orbits));

  /*
  print_clust(
    diff_trans_orbits.begin(),
    diff_trans_orbits.end(),
    std::cout,
    PrototypePrinter<Kinetics::DiffusionTransformation>());
  */

  // Make test config
  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();
  Supercell scel {&primclex, Lattice(2 * a, 2 * b, 3 * c)};
  Configuration config(scel);
  config.init_occupation();
  config.init_displacement();
  config.init_deformation();
  config.init_specie_id();
  config.set_occupation({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0});

  // Make test prim_config
  Configuration prim_config = config.primitive().in_canonical_supercell();

  // IntegralCluster tests
  {
    std::vector<IntegralCluster> scel_generators;
    std::vector<Index> orbit_index;
    std::vector<Index> scel_suborbit_size;
    std::vector<IntegralCluster> config_generators;

    test_1(config, prim_config, orbits, scel_generators, orbit_index, scel_suborbit_size);
    test_2(prim_config, orbits, scel_generators, config_generators, orbit_index, scel_suborbit_size);
    test_3(prim_config, orbits, config_generators);
    test_4(config, orbits);
  }

  // DiffusionTransformation tests
  {
    std::vector<Kinetics::DiffusionTransformation> scel_generators;
    std::vector<Index> orbit_index;
    std::vector<Index> scel_suborbit_size;
    std::vector<Kinetics::DiffusionTransformation> config_generators;

    test_1(config, prim_config, diff_trans_orbits, scel_generators, orbit_index, scel_suborbit_size);
    test_2(prim_config, diff_trans_orbits, scel_generators, config_generators, orbit_index, scel_suborbit_size);
    test_3(prim_config, diff_trans_orbits, config_generators);
    test_4(config, diff_trans_orbits);
  }
}
BOOST_AUTO_TEST_SUITE_END()
