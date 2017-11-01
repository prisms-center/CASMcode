#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/database/Selection.hh"


/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/ConfigEnumAllOccupations_impl.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/container/Enumerator_impl.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/app/enum.hh"
#include "casm/database/Database.hh"
#include "casm/kinetics/DiffusionTransformationEnum.hh"
#include "casm/kinetics/DiffTransConfigEnumOccPerturbations.hh"


using namespace CASM;

BOOST_AUTO_TEST_SUITE(Selection_Test)

BOOST_AUTO_TEST_CASE(Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  const Structure &prim(primclex.prim());
  primclex.settings().set_crystallography_tol(1e-5);

  Completer::EnumOption enum_opt;
  enum_opt.desc();

  // -- Generate Supercell & Configuration --

  ScelEnumByProps enum_scel(primclex, ScelEnumProps(1, 5));
  BOOST_CHECK_EQUAL(true, true);

  ConfigEnumAllOccupations::run(primclex, enum_scel.begin(), enum_scel.end());
  BOOST_CHECK_EQUAL(true, true);


  // Test Supercell Selection
  {
    auto &dict = primclex.settings().query_handler<Supercell>().dict();
    BOOST_CHECK_EQUAL(primclex.generic_db<Supercell>().size(), 13);
    DB::Selection<Supercell> selection(primclex);
    BOOST_CHECK_EQUAL(selection.size(), 13);
    BOOST_CHECK_EQUAL(selection.selected_size(), 0);
  }

  // Test Configuration Selection
  {
    BOOST_CHECK_EQUAL(primclex.generic_db<Configuration>().size(), 126);

    DB::Selection<Configuration> selection(primclex);
    BOOST_CHECK_EQUAL(selection.size(), 126);
    BOOST_CHECK_EQUAL(selection.selected_size(), 0);

    auto &dict = primclex.settings().query_handler<Configuration>().dict();
    selection.set(dict, "lt(scel_size,3)");
    BOOST_CHECK_EQUAL(selection.selected_size(), 9);
  }

  // Test PrimPeriodicDiffTransOrbit Selection
  {
    // Generate
    fs::path difftrans_path = "tests/unit/kinetics/FCCTernary_diff_trans_0.json";
    jsonParser diff_trans_json {difftrans_path};
    Kinetics::DiffusionTransformationEnum::run(primclex, diff_trans_json, enum_opt);
    BOOST_CHECK_EQUAL(true, true);

    // Test
    auto &dict = primclex.settings().query_handler<Kinetics::PrimPeriodicDiffTransOrbit>().dict();
    BOOST_CHECK_EQUAL(primclex.generic_db<Kinetics::PrimPeriodicDiffTransOrbit>().size(), 28);
    DB::Selection<Kinetics::PrimPeriodicDiffTransOrbit> selection(primclex);
    BOOST_CHECK_EQUAL(selection.size(), 28);
    BOOST_CHECK_EQUAL(selection.selected_size(), 0);
  }

  // Test DiffTransConfiguration Selection
  {
    // Generate
    fs::path diffperturb_path = "tests/unit/kinetics/FCCTernary_diff_perturb_0.json";
    jsonParser diff_perturb_json {diffperturb_path};
    Kinetics::DiffTransConfigEnumOccPerturbations::run(primclex, diff_perturb_json, enum_opt);
    BOOST_CHECK_EQUAL(true, true);

    // Test (quantity 1856 not checked for accuracy)
    Kinetics::DiffTransConfigEnumOccPerturbations::run(primclex, diff_perturb_json, enum_opt);
    auto &dict = primclex.settings().query_handler<Kinetics::DiffTransConfiguration>().dict();
    BOOST_CHECK_EQUAL(primclex.generic_db<Kinetics::DiffTransConfiguration>().size(), 1856);
    DB::Selection<Kinetics::DiffTransConfiguration> selection(primclex);
    BOOST_CHECK_EQUAL(selection.size(), 1856);
    BOOST_CHECK_EQUAL(selection.selected_size(), 0);
  }
}

BOOST_AUTO_TEST_SUITE_END()
