#include "gtest/gtest.h"
#include "autotools.hh"

/// What is being tested:
#include "casm/clusterography/ClusterOrbits.hh"

/// What is being used to test it:
#include "casm/clex/PrimClex.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/AppIO_impl.hh"
#include "Common.hh"
#include "ZrOProj.hh"
#include "casm/completer/Handlers.hh"
#include "casm/app/casm_functions.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffusionTransformationEnum.hh"

#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"
#include "casm/symmetry/InvariantSubgroup_impl.hh"

using namespace CASM;
using namespace test;

TEST(LocalClusterExpansionTest, Test0) {

  test::ZrOProj proj;
  proj.check_init();
  proj.check_composition();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);

  fs::path bspecs_path = autotools::abs_srcdir() + "/tests/unit/kinetics/ZrO_bspecs_0.json";
  jsonParser bspecs {bspecs_path};


  std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
  make_prim_periodic_orbits(
    primclex.shared_prim(),
    bspecs,
    alloy_sites_filter,
    primclex.crystallography_tol(),
    std::back_inserter(orbits),
    primclex.log());

  std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
  Kinetics::make_prim_periodic_diff_trans_orbits(
    orbits.begin() + 2,
    orbits.begin() + 3,
    primclex.crystallography_tol(),
    std::back_inserter(diff_trans_orbits),
    &primclex);
  Kinetics::DiffusionTransformation trans = diff_trans_orbits[0].prototype();

  fs::path local_bspecs_path = autotools::abs_srcdir() + "/tests/unit/kinetics/ZrO_local_bspecs_1.json";
  jsonParser local_bspecs {local_bspecs_path};

  // in ZrO we are looking at the hop from octahedral interstitial (O-site) to nearest octahedral interstitial (O-site):

  // THE ORDERING OF THESE ORBITS MAY CHANGE IF THE ENUMERATION PROCESS IS SORTED BY dist_to_path()
  // empty cluster has 1 option
  // point cluster has 6 options given the cut off radius of 6 angstroms and ZrO unitcell params
  //    2 along c axis (1 2.59 Angs away 1 5.18 Angs away) multiplicity 2 #1-2
  //    2 that make hexagonal prisms with their equivalents (1 3.23 Angs away, 1 4.145 Angs away) multiplicity 12 #3-4
  //    2 that make triangular prisms with their equivalents (both 5.60 Angs away rotated 60 degrees from each other) multiplicity 6 #5-6
  // 2-cluster has 26 options given the cut off radius of 6 angstroms and ZrO unitcell params
  //  a #1-#2 cluster on one side of the trans Mult: 2  Length: 2.59 #7
  //    #1-#4 cluster on one side of the trans Mult: 12 Length: 3.23 #12
  //    #1-#3 cluster on one side of the trans Mult: 12 Length: 4.145 #19
  //    #2-#4 cluster on one side of the trans Mult: 12 Length: 4.145 #20
  //    #3-#3 (hop translation) Mult: 6 Length: 2.59 #8
  //    #3-#3 perpendicular to hop same c value base towards #5 Mult: 6 Length: 3.23 #13
  //    #3-#3 perpendicular to hop same c value base towards #6 Mult: 6 Length: 3.23 #15
  //    #3-#3 perpendicular to hop same c value Mult: 12 Length: 5.60 #30
  //    #3-#3 slanted to hop diff c value (left handed) Mult: 6 Length: 4.145 #21
  //    #3-#3 slanted to hop diff c value (right handed) Mult: 6 Length: 4.145 #23
  //    #3-#4 shifted 0.5 c from hop translation Mult:12 Length 2.59 #9
  //    #3-#4 slanted to hop (left handed slant) Mult: 12 Length: 4.145 #22
  //    #3-#4 slanted to hop (right handed slant) Mult: 12 Length: 4.145 #24
  //    #3-#4 parallel to hop Mult: 12 Length: 5.18 #29
  //    #3-#5 perpendicular to hop same c value Mult: 12 Length: 3.23 #17
  //    #3-#5 slanted to hop  Mult: 12 Length: 4.145 #25
  //    #3-#6 slanted to hop Mult: 12 Length: 4.145 #26
  //    #3-#6 perpendicular to hop same c value Mult: 12 Length: 3.23 #18
  //    #4-#4 perpendicular to hop same c value base towards #5 Mult: 6 Length: 3.23 #14
  //    #4-#4 perpendicular to hop same c value base towards #6 Mult: 6 Length: 3.23 #16
  //    #4-#4 perpendicular to hop same c value Mult: 12 Length: 5.60 #31
  //    #4-#5 slanted to hop  Mult: 12 Length: 4.145 #27
  //    #4-#6 slanted to hop Mult: 12 Length: 4.145 #28
  //    #5-#5 (hop translation) Mult: 3 Length 2.59 #10
  //    #5-#6 perpendicular to hop Mult: 12 Length 5.60 #32
  //    #6-#6 (hop translation) Mult: 3 Length 2.59 #11
  // 3-cluster have 4 options given the cut off radius of 5 angstroms a max cluster size of 4 angstroms and ZrO unitcell params
  //    #1-#4-#4 base towards #5 Mult: 6 MinLength 3.23 MaxLength 3.23 #33
  //    #1-#4-#4 base towards #6 Mult: 6 MinLength 3.23 MaxLength 3.23 #34
  //    #3-#3-#5 Mult: 6 MinLength 3.23 MaxLength 3.23 #35
  //    #3-#3-#6 Mult: 6 MinLength 3.23 MaxLength 3.23 #36

  std::map<int, int> count_check = {{0, 1}, {1, 6}, {2, 26}, {3, 4}};



  // these have been manually checked
  std::vector<int> mult;
  std::vector<int> expected_mult = {1, 2, 2, 12, 12, 6, 6, 2, 6, 12, 3, 3, 12, 6, 6, 6, 6, 12, 12, 12, 12, 6, 12, 6, 12,
                                    12, 12, 12, 12, 12, 12, 12, 12, 6, 6, 6, 6
                                   };
  auto expected_mult_it = expected_mult.begin();

  std::vector<LocalIntegralClusterOrbit> local_orbits;
  SymGroup generating_grp {
    trans.invariant_subgroup(
      primclex.prim().factor_group(),
      PrimPeriodicDiffTransSymCompare(primclex.shared_prim(), primclex.crystallography_tol()))};
  LocalSymCompare<IntegralCluster> sym_compare(primclex.shared_prim(), primclex.crystallography_tol());

  make_local_orbits(
    trans,
    generating_grp,
    sym_compare,
    local_bspecs,
    alloy_sites_filter,
    primclex.crystallography_tol(),
    std::back_inserter(local_orbits),
    primclex.log());


  // check the size of the local orbits
  for(const auto &orb : local_orbits) {
    mult.push_back(orb.size());
    if(expected_mult_it != expected_mult.end()) {
      EXPECT_EQ(orb.size(), *expected_mult_it);
      ++expected_mult_it;
      --count_check[orb.prototype().size()];
    }
    else {
      EXPECT_EQ(orb.size(), 0);
    }
  }

  for(const auto &pair : count_check) {
    EXPECT_EQ(pair.second, 0);
  }


  //std::cout << trans << std::endl;
  //PrototypePrinter<IntegralCluster> printer;
  //print_clust(local_orbits.begin(), local_orbits.end(), std::cout, printer);



}
