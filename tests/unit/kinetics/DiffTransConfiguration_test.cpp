#include "gtest/gtest.h"
#include "autotools.hh"

/// What is being tested:
#include "casm/kinetics/DiffTransConfiguration_impl.hh"

/// What is being used to test it:
#include "casm/clex/PrimClex.hh"
#include "casm/app/AppIO.hh"
#include "Common.hh"
#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"
#include "casm/clusterography/ClusterOrbits.hh"

using namespace CASM;
using namespace test;

TEST(DiffTransConfigurationTest, Test0) {

  test::ZrOProj proj;
  proj.check_init();
  proj.check_composition();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);

  fs::path bspecs_path = autotools::abs_srcdir() + "/tests/unit/kinetics/ZrO_bspecs_0.json";
  jsonParser bspecs {bspecs_path};

  std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
  make_prim_periodic_orbits(
    primclex.prim(),
    bspecs,
    alloy_sites_filter,
    primclex.crystallography_tol(),
    std::back_inserter(orbits),
    primclex.log());

  //print_clust(orbits.begin(), orbits.end(), primclex.log(), ProtoSitesPrinter());
  std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
  Kinetics::make_prim_periodic_diff_trans_orbits(
    orbits.begin() + 4,
    orbits.begin() + 7,
    primclex.crystallography_tol(),
    std::back_inserter(diff_trans_orbits),
    &primclex);
  Kinetics::DiffusionTransformation trans = diff_trans_orbits[0].prototype();
  Kinetics::DiffusionTransformation trans2 = diff_trans_orbits[2].prototype();

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  Supercell scel {&primclex, Lattice(2 * a, 2 * b, 2 * c)};

  Configuration config = Configuration::zeros(scel);
  //hardcoded occupation for trans to occur is there a way to do this generally?
  config.set_occupation(std::vector<int>({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}));

  Configuration config2 = Configuration::zeros(scel);
  //hardcoded occupation for trans to occur is there a way to do this generally?
  config2.set_occupation(std::vector<int>({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0}));

  Configuration config3 = Configuration::zeros(scel);
  config3.set_occupation(std::vector<int>({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}));

  //test make attachable
  Configuration result = make_attachable(trans, config3);
  EXPECT_EQ(config3 == result, 0);
  for(auto traj : trans.species_traj()) {
    Index l = result.supercell().linear_index(traj.from.uccoord);
    EXPECT_EQ(result.occ(l), traj.from.occ);
  }

  //test Constructor/field accessors
  Kinetics::DiffTransConfiguration dtc(make_attachable(trans, config), trans);
  EXPECT_EQ(dtc.from_config(), make_attachable(trans, config));
  EXPECT_EQ(dtc.diff_trans(), trans);
  Configuration tmp {config};
  tmp = make_attachable(trans, tmp);
  tmp = dtc.diff_trans().apply_to(tmp);
  EXPECT_EQ(dtc.to_config(), tmp);

  //check comparison
  Kinetics::DiffTransConfiguration dtc2(make_attachable(trans2, config2), trans2);
  Kinetics::DiffTransConfiguration dtc3(make_attachable(trans, config2), trans);

  //config > config2 but trans < trans2
  //comparing transformation takes priority
  EXPECT_EQ(dtc < dtc2, trans < trans2);
  //If the above test fails the preference of priority may have changed
  EXPECT_EQ(dtc < dtc3, config < config2);

  //check apply sym
  PermuteIterator it = config.supercell().sym_info().permute_begin();
  EXPECT_EQ(copy_apply(it, dtc) == dtc, 1);

  it = it.begin_next_fg_op();
  Configuration new_config = copy_apply(it, make_attachable(trans, config));
  Kinetics::ScelPeriodicDiffTransSymCompare symcompare(config.supercell().prim_grid(),
                                                       config.supercell().crystallography_tol());
  Kinetics::DiffusionTransformation new_trans =
    symcompare.prepare(copy_apply(it.sym_op(), trans));

  //std::cout << newdtc << std::endl;

  EXPECT_EQ(copy_apply(it, dtc).from_config() == new_config, 1);
  EXPECT_EQ(copy_apply(it, dtc).diff_trans() == new_trans, 1);

  //check sorting
  EXPECT_EQ(dtc.is_sorted(), dtc.from_config() < dtc.to_config());
  EXPECT_EQ(dtc.is_sorted(), dtc == dtc.sorted());


  //check canonical form
  EXPECT_EQ(dtc.is_canonical(), 0);

  EXPECT_EQ(!dtc.is_canonical(), dtc < dtc.canonical_form());
  EXPECT_EQ(dtc.is_canonical(), dtc == dtc.canonical_form());
  EXPECT_EQ(1, dtc.canonical_form().is_canonical());
  EXPECT_EQ(1, copy_apply(dtc.to_canonical(), dtc) == dtc.canonical_form());
  EXPECT_EQ(1, copy_apply(dtc.from_canonical(), dtc.canonical_form()) == dtc);
  jsonParser dtcjson;
  dtcjson.put_obj();
  dtc.to_json(dtcjson);
  Kinetics::DiffTransConfiguration loaded_dtc(dtc.from_config().supercell(), dtcjson);
  EXPECT_EQ(dtc == loaded_dtc, 1);

}

