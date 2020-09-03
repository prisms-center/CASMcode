#include "gtest/gtest.h"

#include "casm/casm_io/Log.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"

#include "crystallography/TestStructures.hh" // for test::ZrO_prim

// An "orbit" is the set of equivalent elements that can be generated under application of
//   symmetry operations in a group. In CASM, an Orbit instance is constructed from:
//   - a "generating element", one element which will be transformed by operations in a symmetry
//     group generating all the equivalent elements in the orbit
//   - a "generating group", the symmetry group to be applied
//   - a "SymCompare" functor, which customizes the orbit constructing algorithm to be able to
//     apply symmetry transformation to elements and efficiently compare them:
//     - SymCompare::copy_apply: apply symmetry to generate new elements
//     - SymCompare::spatial_prepare: transform an element to a standard spatial position
//       (i.e. translate into the unit cell) for efficient comparison
//     - SymCompare::representation_prepare: transform an element to a standard representation for
//       efficient comparison (i.e. sort cluster elements)
//     - SymCompare::prepare: prepare an element for comparison via representation_prepare,
//       followed by spatial_prepare
//     - SymCompare::compare: order prepared elements
//     - SymCompare::make_invariants: construct orbit invariants from one element in the orbit,
//       for efficient comparison (i.e. cluster size, site-to-site distances, etc.)
//     - SymCompare::compare_invariants: order orbits by invariants
//     - SymCompare::inter_orbit_compare: order orbits, breaking invariants_compare ties
//     - SymCompare::canonical_transform: access the transform that took an element from
//       the unprepared to prepared state
//     - SymCompare::spatial_transform: access the spatial transform that was used during most
//       recent spatial preparation of an element
//
// CASM provides SymCompare functors for several standard cluster orbit uses cases:
// - PrimPeriodicSymCompare: For comparing elements with the translation symmetry of the "prim"
//   structure's lattice. This is the appropriate method for generating the cluster expansion
//   basis for crystal properties.
// - ScelPeriodicSymCompare: For comparing elements with the translation symmetry of a superlattice
//   of the "prim" structure's lattice. This is the method that is analogous to
//   PrimPeriodicSymCompare, but for a non-primitive structure.
// - AperiodicSymCompare (or alias LocalSymCompare): For comparing elements without translation
//   symmetry. This is the appropriate method for generating the cluster expansion basis for
//   properties local to a site or cluster of sites in a crystal.
// - WithScelSymCompare: For comparing elements by placing all sites "within" a supercell and
//   then comparing. When using this method there are a finite number of orbits that may possibly
//   be generated rather than an infinite number. For example, sites in a pair cluster cannot get
//   separated farther than the half the maximum periodic image distance because they are returned
//   "within" the supercell. This is the appropriate method for generating the unique clusters
//   in a supercell with periodic boundary conditions.
//
// After the Orbit is constructed it has:
// - an "elements" vector of all elements in the orbit
// - a "prototype" element, the first element in the elements vector
// - an "equivalence map", a table of SymOp such that equivalence_map()[i][j], for all j, acts
//   on the prototype element to generate element(i). More specifically:
//     equivalence_map()[i][j] = t * g, where
//     g is a generating group element CASM::SymOp
//     t is the "spatial_transform" CASM::SymOp such that
//       sym_compare.representation_prepare(element(i)) ==
//       sym_compare.copy_apply(t, sym_compare.representation_prepare(sym_compare.copy_apply(g, prototype)))
//   The first row of the equivalence map is the subgroup of the generating group that leaves the
//   prototype invariant.

// Function asserting attributes true of all orbits
template<typename OrbitType>
void expect_for_any_orbit(
  OrbitType const &orbit,
  std::shared_ptr<CASM::Structure const> const &shared_prim,
  std::string test_case) {

  auto const &sym_compare = orbit.sym_compare();

  // the orbit size is the orbit's element vector size
  EXPECT_EQ(orbit.size(), orbit.elements().size()) << test_case + ": orbit size error";

  // the prototype is "prepared"
  EXPECT_TRUE(sym_compare.equal(sym_compare.prepare(orbit.prototype()), orbit.prototype())) <<
      test_case + ": prototype preparation error";

  // the "prototype" is the first elment in the orbit's element vector
  for(int i = 0; i < orbit.size(); ++i) {
    EXPECT_EQ(sym_compare.equal(orbit.prototype(), orbit.element(i)), i == 0) <<
                                                                              test_case + ": prototype location error";
  }

  // Check the equivalence_map

  // # of rows == # elements
  EXPECT_EQ(orbit.size(), orbit.equivalence_map().size()) <<
                                                          test_case + ": number of equivalence_map rows error";

  for(int i = 0; i < orbit.equivalence_map().size(); ++i) {

    // number of columns == generating group size / number of elements
    EXPECT_EQ(orbit.equivalence_map()[i].size(), orbit.generating_group().size() / orbit.size()) <<
        test_case + ": number of equivalence_map columns error";

    for(int j = 0; j < orbit.equivalence_map()[i].size(); ++j) {
      // apply symop from the equivalence_map to prototype
      auto cluster = CASM::sym::copy_apply(
                       orbit.equivalence_map()[i][j], orbit.prototype(), *shared_prim);

      // for first column only, applying equivalence_map symop to prototype yields element(i) directly
      if(j == 0) {
        EXPECT_TRUE(sym_compare.equal(orbit.element(i), cluster)) <<
                                                                  test_case + ": not equivalent after applying first equivalence map symop";
      }

      // for other columns, applying equivalence_map symop may result in a different representation
      auto rep_prepared_element_i = sym_compare.representation_prepare(orbit.element(i));
      auto rep_prepared_cluster = sym_compare.representation_prepare(cluster);
      EXPECT_TRUE(sym_compare.equal(rep_prepared_element_i, rep_prepared_cluster)) <<
                                                                                   test_case + ": not equivalent after representation prepare";

      // result will already be "spatially prepared"
      //   (representation prepare is always done before spatial prepare)
      auto prepared_cluster = sym_compare.spatial_prepare(rep_prepared_cluster);
      EXPECT_TRUE(sym_compare.equal(prepared_cluster, rep_prepared_cluster)) <<
                                                                             test_case + ": not equivalent after spatial prepare";

    }
  }
}

TEST(ExampleClusterographyClusterOrbits, OrbitConstructorTernaryFCC) {

  // Example constructing a cluster orbit:

  // Construct the prim structure: FCC, with [A, B, C] allowed occupants
  // - Lattice:
  //   - a: [0.0, 2.0, 2.0]
  //   - b: [2.0, 0.0, 2.0]
  //   - c: [2.0, 2.0, 0.0]
  // - Basis:
  //   - 0: (0., 0., 0.), [A, B, C]
  auto shared_prim = std::make_shared<CASM::Structure const>(test::FCC_ternary_prim());

  // Tolerance value to use
  double xtal_tol = shared_prim->lattice().tol();

  {
    // Construct a pair cluster of two [A, B, C] sites, the nearest neighbors
    //   There are 6 equivalent clusters of this type
    CASM::IntegralCluster cluster {*shared_prim};
    cluster.elements().emplace_back(0, 0, 0, 0);
    cluster.elements().emplace_back(0, 1, 0, 0);

    // Construct an orbit using PrimPeriodicSymCompare and the prim factor group
    typedef CASM::PrimPeriodicSymCompare<CASM::IntegralCluster> sym_compare_type;
    typedef CASM::Orbit<sym_compare_type> orbit_type;
    orbit_type orbit {cluster, shared_prim->factor_group(), sym_compare_type {shared_prim, xtal_tol}};

    // --- The following are true for this example ---
    EXPECT_EQ(shared_prim->factor_group().size(), 48);
    EXPECT_EQ(orbit.size(), 6);

    // --- The following function asserts things that are true for any orbit ---
    expect_for_any_orbit(orbit, shared_prim, "FCC nearest neighbor pair");

  }

  {
    // Construct a pair cluster of two [A, B, C] sites, the conventional cubic cell edges
    //   There are 3 equivalent clusters or this type, corresponding to x, y, z axes
    CASM::IntegralCluster cluster {*shared_prim};
    cluster.elements().emplace_back(0, 0, 0, 0);
    cluster.elements().emplace_back(0, 1, 1, -1);

    // Construct an orbit using PrimPeriodicSymCompare and the prim factor group
    typedef CASM::PrimPeriodicSymCompare<CASM::IntegralCluster> sym_compare_type;
    typedef CASM::Orbit<sym_compare_type> orbit_type;
    orbit_type orbit {cluster, shared_prim->factor_group(), sym_compare_type {shared_prim, xtal_tol}};

    // --- The following are true for this example ---
    EXPECT_EQ(shared_prim->factor_group().size(), 48);
    EXPECT_EQ(orbit.size(), 3);

    // --- The following function asserts things that are true for any orbit ---
    expect_for_any_orbit(orbit, shared_prim, "FCC conventional cell edge pair");

  }

}

TEST(ExampleClusterographyClusterOrbits, OrbitConstructorZrO) {

  // Example constructing a cluster orbit:

  // Construct the prim structure: HCP Zr, with [Va, O] octahedral interstitial sites
  // - Lattice:
  //   - a: [3.234, -1.617, 0.000]
  //   - b: [0.000,  2.801, 0.000]
  //   - c: [0.000,  0.000, 5.169]
  // - Basis:
  //   - 0: (0., 0., 0.), [Zr]
  //   - 1: (2./3., 1./3., 1./2.), [Zr]
  //   - 2: (1./3., 2./3., 1./4.), [Va, O]
  //   - 3: (1./3., 2./3., 3./4.), [Va, O]
  auto shared_prim = std::make_shared<CASM::Structure const>(test::ZrO_prim());

  // Tolerance value to use
  double xtal_tol = shared_prim->lattice().tol();

  {
    // Construct a pair cluster of two [Va, O] sites, the closest neighbors along the c-axis
    //   There are 2 equivalent clusters of this type
    CASM::IntegralCluster cluster {*shared_prim};
    cluster.elements().emplace_back(2, 0, 0, 0);
    cluster.elements().emplace_back(3, 0, 0, 0);

    // Construct an orbit using PrimPeriodicSymCompare and the prim factor group
    typedef CASM::PrimPeriodicSymCompare<CASM::IntegralCluster> sym_compare_type;
    typedef CASM::Orbit<sym_compare_type> orbit_type;
    orbit_type orbit {cluster, shared_prim->factor_group(), sym_compare_type {shared_prim, xtal_tol}};

    // --- The following are true for this example ---
    EXPECT_EQ(shared_prim->factor_group().size(), 24);
    EXPECT_EQ(orbit.size(), 2);

    // --- The following function asserts things that are true for any orbit ---
    expect_for_any_orbit(orbit, shared_prim, "ZrO c-axis pair");

  }

  {
    // Construct a pair cluster of two [Va, O] sites, the closest neighbors along the a-axis
    //   There are 6 equivalent clusters of this type, 3 each involving basis site 2 and basis site 3
    CASM::IntegralCluster cluster {*shared_prim};
    cluster.elements().emplace_back(2, 0, 0, 0);
    cluster.elements().emplace_back(2, 1, 0, 0);

    // Construct an orbit using PrimPeriodicSymCompare and the prim factor group
    typedef CASM::PrimPeriodicSymCompare<CASM::IntegralCluster> sym_compare_type;
    typedef CASM::Orbit<sym_compare_type> orbit_type;
    orbit_type orbit {cluster, shared_prim->factor_group(), sym_compare_type {shared_prim, xtal_tol}};

    // --- The following are true for this example ---
    EXPECT_EQ(shared_prim->factor_group().size(), 24);
    EXPECT_EQ(orbit.size(), 6);

    int basis_site_2_count = 0;
    int basis_site_3_count = 0;
    for(auto const &cluster : orbit.elements()) {
      EXPECT_EQ(cluster[0].sublattice(), cluster[1].sublattice());
      if(cluster[0].sublattice() == 2) {
        basis_site_2_count++;
      }
      if(cluster[0].sublattice() == 3) {
        basis_site_3_count++;
      }
    }
    EXPECT_EQ(basis_site_2_count, 3);
    EXPECT_EQ(basis_site_3_count, 3);

    // --- The following function asserts things that are true for any orbit ---
    expect_for_any_orbit(orbit, shared_prim, "ZrO a-axis pair");

  }

}

TEST(ExampleClusterographyClusterOrbits, OrbitGenerators) {
  // Multiple orbits may be generated using the OrbitGenerators class to identify a set of unique
  //   generating elements and then use those to construct unique orbits.

  // Construct a ZrO prim (same as above)
  auto shared_prim = std::make_shared<CASM::Structure const>(test::ZrO_prim());

  // Tolerance value to use
  double xtal_tol = shared_prim->lattice().tol();

  typedef CASM::PrimPeriodicSymCompare<CASM::IntegralCluster> sym_compare_type;
  typedef CASM::Orbit<sym_compare_type> orbit_type;

  // Construct a sym_compare functor
  sym_compare_type sym_compare {shared_prim, xtal_tol};

  // Construct the OrbitGenerators instance with the generating group and sym_compare
  CASM::OrbitGenerators<orbit_type> generators {shared_prim->factor_group(), sym_compare};

  // Generate orbits for all single site clusters in the prim unit cell:

  // First, insert all single site clusters into generators, which will keep unique clusters
  CASM::IntegralCluster test_cluster {*shared_prim};
  for(int b = 0; b < shared_prim->basis().size(); ++b) {
    test_cluster.elements().clear();
    test_cluster.elements().emplace_back(b, 0, 0, 0);

    // The `insert` member will use the sym_compare object to transform test_cluster into its
    //   "prepared" form and insert it into a set of unique generating elements
    // There also exists an `insert_canonical` member which can be used if the cluster to be
    //   inserted has already been "prepared" externally.
    generators.insert(test_cluster);
  }

  // there are two unique clusters in the ZrO prim: a [Zr] basis site, and a [Va, O] basis site
  EXPECT_EQ(generators.elements.size(), 2);


  // Second, make orbits from the unique clusters

  // Use the generating elements in generators.elements to construct orbits
  std::vector<orbit_type> orbits;
  generators.make_orbits(std::back_inserter(orbits));

  // From the 2 generating clusters there are two orbits, each with 2 elements
  EXPECT_EQ(generators.elements.size(), 2);
  EXPECT_EQ(orbits[0].size(), 2);
  EXPECT_EQ(orbits[1].size(), 2);
}

TEST(ExampleClusterographyClusterOrbits, AsymmetricUnitOrbits) {

  // The asymmetric unit is the set of symmetrically unique sites in the unit cell of the prim.
  //
  // The asymmetric unit is important for structural analysis and for constructing the cluster
  //   expansion basis. Site basis functions must be constructed for each prototype site in the
  //   asymmetric unit, and then each symmetrically equivalent site must have symmetrically
  //   equivalent site basis functions.

  // Construct a ZrO prim (same as above)
  auto shared_prim = std::make_shared<CASM::Structure const>(test::ZrO_prim());

  // Tolerance value to use
  double xtal_tol = shared_prim->lattice().tol();

  // Container that will hold resulting orbits
  //   (PrimPeriodicIntegralClusterOrbit is a typedef for Orbit<PrimPeriodicSymCompare<IntegralCluster>>)
  std::vector<CASM::PrimPeriodicIntegralClusterOrbit> asymmetric_unit_orbits;

  // A function with signature `bool (xtal::Site)` which returns True for sites that should be
  //   included in the generated orbits, and False for sites that will be excluded. Options:
  //   - all_sites_filter: include all sites
  //   - alloy_sites_filter: include sites with >1 allowed occupant DoF
  //   - dof_sites_filter(const std::vector<DoFKey> &dofs = {}):
  //     - If dofs is empty (default), include sites that have any continuous DoF or >1 allowed
  //       occupant DoF
  //     - If dofs is not empty, include sites that have any of the DoF types included. Use "occ"
  //       to include sites with >1 allowed occupant DoF.
  //
  // In the examples below this will be varied.
  CASM::SiteFilterFunction site_filter;

  // Orbit generation can be slow depending on the generating criteria, so orbit generating functions
  //   in CASM typically take a `status` argument of type `std::ostream&` where status info can be
  //   printed while the orbit generation method is running. Asymmetric units are generally small
  //   and quick, so here we pass the CASM::null_log() and nothing will be printed.
  std::ostream &status = CASM::null_log();

  {
    // Include all sites:
    // - Expect two orbits: [Zr] sites and [Va,O] sites
    site_filter = CASM::all_sites_filter;

    asymmetric_unit_orbits.clear();
    auto inserter = std::back_inserter(asymmetric_unit_orbits);
    make_prim_periodic_asymmetric_unit(shared_prim, site_filter, xtal_tol, inserter, status);
    EXPECT_EQ(asymmetric_unit_orbits.size(), 2);
  }

  {
    // Include alloying sites:
    // - Expect one orbit: [Va,O] sites
    site_filter = CASM::alloy_sites_filter;

    asymmetric_unit_orbits.clear();
    auto inserter = std::back_inserter(asymmetric_unit_orbits);
    make_prim_periodic_asymmetric_unit(shared_prim, site_filter, xtal_tol, inserter, status);
    EXPECT_EQ(asymmetric_unit_orbits.size(), 1);
  }

  {
    // Include "occ" sites:
    // - Expect one orbit: [Va,O] sites
    site_filter = CASM::dof_sites_filter({"occ"});

    asymmetric_unit_orbits.clear();
    auto inserter = std::back_inserter(asymmetric_unit_orbits);
    make_prim_periodic_asymmetric_unit(shared_prim, site_filter, xtal_tol, inserter, status);
    EXPECT_EQ(asymmetric_unit_orbits.size(), 1);
  }
}
