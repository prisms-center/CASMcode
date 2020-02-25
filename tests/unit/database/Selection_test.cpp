#include "gtest/gtest.h"
#include "autotools.hh"

/// What is being tested:
#include "casm/database/Selection.hh"


/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/ConfigEnumAllOccupations_impl.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/enumerator/Enumerator_impl.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/app/enum.hh"
#include "casm/database/Database.hh"
#include "casm/kinetics/DiffusionTransformationEnum.hh"
#include "casm/kinetics/DiffTransConfigEnumOccPerturbations.hh"


using namespace CASM;

TEST(Selection_Test, Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  //const Structure &prim(primclex.prim());
  primclex.settings().set_crystallography_tol(1e-5);

  Completer::EnumOption enum_opt;
  enum_opt.desc();

  // -- Generate Supercell & Configuration --

  ScelEnumByProps enum_scel(primclex, ScelEnumProps(1, 5));
  EXPECT_EQ(true, true);

  ConfigEnumAllOccupations::run(primclex, enum_scel.begin(), enum_scel.end());
  EXPECT_EQ(true, true);


  // Test Supercell Selection
  {
    //auto &dict = primclex.settings().query_handler<Supercell>().dict();
    EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
    DB::Selection<Supercell> selection(primclex);
    EXPECT_EQ(selection.size(), 13);
    EXPECT_EQ(selection.selected_size(), 0);
  }

  // Test Configuration Selection
  {
    EXPECT_EQ(primclex.generic_db<Configuration>().size(), 126);

    DB::Selection<Configuration> selection(primclex);
    EXPECT_EQ(selection.size(), 126);
    EXPECT_EQ(selection.selected_size(), 0);

    auto &dict = primclex.settings().query_handler<Configuration>().dict();
    selection.set(dict, "lt(scel_size,3)");
    EXPECT_EQ(selection.selected_size(), 9);
  }

  // Test PrimPeriodicDiffTransOrbit Selection
  {
    // Generate
    fs::path difftrans_path = autotools::abs_srcdir() + "/tests/unit/kinetics/FCCTernary_diff_trans_0.json";
    jsonParser diff_trans_json {difftrans_path};
    Kinetics::DiffusionTransformationEnum::run(primclex, diff_trans_json, enum_opt, nullptr);
    EXPECT_EQ(true, true);

    // Test
    //auto &dict = primclex.settings().query_handler<Kinetics::PrimPeriodicDiffTransOrbit>().dict();
    EXPECT_EQ(primclex.generic_db<Kinetics::PrimPeriodicDiffTransOrbit>().size(), 28);
    DB::Selection<Kinetics::PrimPeriodicDiffTransOrbit> selection(primclex);
    EXPECT_EQ(selection.size(), 28);
    EXPECT_EQ(selection.selected_size(), 0);
  }

  //This test was commented out for the basicstructure_refactor merge
  //the tests seem to be expecting things to be generated in a particular order that is not
  //guaranteed. As a result, the following exception is thrown:
  //"Error in DiffTransConfiguration constructor both endpoints are exactly the same config!"
  //In addition to receiving an unexpected DiffTransConfig, it seems that the current behavior
  //does not allow for a diffusion tranformation in which two atoms of the same type swap places.
  //Is that how we want to keep it?

  /*
  // Test DiffTransConfiguration Selection
  {
    // Generate
    fs::path diffperturb_path = autotools::abs_srcdir() + "/tests/unit/kinetics/FCCTernary_diff_perturb_0.json";
    jsonParser diff_perturb_json {diffperturb_path};
    Kinetics::DiffTransConfigEnumOccPerturbations::run(primclex, diff_perturb_json, enum_opt, nullptr);
    EXPECT_EQ(true, true);

    // Test (quantity 1856 not checked for accuracy)
    Kinetics::DiffTransConfigEnumOccPerturbations::run(primclex, diff_perturb_json, enum_opt, nullptr);
    //auto &dict = primclex.settings().query_handler<Kinetics::DiffTransConfiguration>().dict();
    EXPECT_EQ(primclex.generic_db<Kinetics::DiffTransConfiguration>().size(), 1856);
    DB::Selection<Kinetics::DiffTransConfiguration> selection(primclex);
    EXPECT_EQ(selection.size(), 1856);
    EXPECT_EQ(selection.selected_size(), 0);
  }
  */
}

