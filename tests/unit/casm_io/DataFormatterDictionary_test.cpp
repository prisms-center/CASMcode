#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/casm_io/DataFormatter.hh"

#include <sstream>
#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/ConfigEnumAllOccupations_impl.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/database/Selection.hh"
#include "casm/database/Selected.hh"
#include "casm/container/Enumerator_impl.hh"
#include "casm/app/enum.hh"
#include "casm/kinetics/DiffusionTransformationEnum.hh"
#include "casm/kinetics/DiffTransConfigEnumOccPerturbations.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(DataFormatterDictionaryTest)

BOOST_AUTO_TEST_CASE(Test1) {

  BOOST_CHECK_EQUAL(true, true);
  test::FCCTernaryProj proj;
  proj.check_init();
  proj.check_composition();

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

  // -- Check Supercell --

  {
    BOOST_CHECK_EQUAL(primclex.db<Supercell>().size(), 13);
    DB::Selection<Supercell> selection(primclex, "ALL");

    auto &qh = primclex.settings().query_handler<Supercell>();
    qh.set_selected(selection);
    auto &dict = qh.dict();

    std::vector<std::string> all_col = {"name", "selected", "scel_size"};
    auto formatter = dict.parse(all_col);

    std::stringstream ss;
    ss << formatter(selection.all().begin(), selection.all().end()) << std::endl;
    BOOST_CHECK_EQUAL(true, true);
  }

  // -- Check Configuration --

  {
    BOOST_CHECK_EQUAL(primclex.db<Configuration>().size(), 126);

    DB::Selection<Configuration> selection(primclex, "ALL");

    auto &qh = primclex.settings().query_handler<Configuration>();
    qh.set_selected(selection);
    auto &dict = qh.dict();

    std::vector<std::string> all_col = {"name", "selected", "comp"};
    auto formatter = dict.parse(all_col);

    std::stringstream ss;
    ss << formatter(selection.all().begin(), selection.all().end()) << std::endl;
    BOOST_CHECK_EQUAL(true, true);
  }
  // -- DiffusionTransformation --

  {
    fs::path difftrans_path = "tests/unit/kinetics/diff_trans.json";
    jsonParser diff_trans_json {difftrans_path};
    BOOST_CHECK_EQUAL(true, true);
    Kinetics::DiffusionTransformationEnum::run(primclex, diff_trans_json, enum_opt);
    BOOST_CHECK_EQUAL(true, true);

    BOOST_CHECK_EQUAL(primclex.db<Kinetics::PrimPeriodicDiffTransOrbit>().size(), 28);
    primclex.db<Kinetics::PrimPeriodicDiffTransOrbit>().commit();
    DB::Selection<Kinetics::PrimPeriodicDiffTransOrbit> selection(primclex, "ALL");

    auto &qh = primclex.settings().query_handler<Kinetics::PrimPeriodicDiffTransOrbit>();
    qh.set_selected(selection);
    auto &dict = qh.dict();

    std::vector<std::string> all_col = {"name", "selected"};
    auto formatter = dict.parse(all_col);

    std::stringstream ss;
    ss << formatter(selection.all().begin(), selection.all().end()) << std::endl;
    BOOST_CHECK_EQUAL(true, true);
  }

  // -- DiffTransConfiguration --

  {
    fs::path diffperturb_path = "tests/unit/kinetics/diff_perturb.json";
    jsonParser diff_perturb_json {diffperturb_path};
    Kinetics::DiffTransConfigEnumOccPerturbations::run(primclex, diff_perturb_json, enum_opt);

    BOOST_CHECK_EQUAL(primclex.db<Kinetics::DiffTransConfiguration>().size(), 2);
    primclex.db<Kinetics::DiffTransConfiguration>().commit();
    DB::Selection<Kinetics::DiffTransConfiguration> selection(primclex, "ALL");

    auto &qh = primclex.settings().query_handler<Kinetics::DiffTransConfiguration>();
    qh.set_selected(selection);
    auto &dict = qh.dict();

    std::vector<std::string> all_col = {"name", "selected"};
    auto formatter = dict.parse(all_col);

    std::stringstream ss;
    ss << formatter(selection.all().begin(), selection.all().end()) << std::endl;
    BOOST_CHECK_EQUAL(true, true);
  }

}

BOOST_AUTO_TEST_SUITE_END()




