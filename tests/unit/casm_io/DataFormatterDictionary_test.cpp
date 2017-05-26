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

using namespace CASM;

BOOST_AUTO_TEST_SUITE(DataFormatterDictionaryTest)

BOOST_AUTO_TEST_CASE(Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();
  proj.check_composition();

  PrimClex primclex(proj.dir, null_log());
  const Structure &prim(primclex.prim());
  primclex.settings().set_crystallography_tol(1e-5);

  fs::path difftrans_path = "tests/unit/kinetics/diff_trans.json";
  jsonParser diff_trans_json {difftrans_path};
  Completer::EnumOption enum_opt;
  enum_opt.desc();

  ScelEnumByProps enum_scel(primclex, ScelEnumProps(1, 5));
  ConfigEnumAllOccupations::run(primclex, enum_scel.begin(), enum_scel.end());
  Kinetics::DiffusionTransformationEnum::run(primclex, diff_trans_json, enum_opt);
  BOOST_CHECK_EQUAL(primclex.db<Configuration>().size(), 126);
  primclex.db<Configuration>().commit();

  DB::Selection<Configuration> selection(primclex, "ALL");

  auto &qh = primclex.settings().query_handler<Configuration>();
  qh.set_selected(selection);
  auto &dict = qh.dict();

  std::vector<std::string> all_col = {"name", "selected", "comp"};
  auto formatter = dict.parse(all_col);

  std::stringstream ss;
  ss << formatter(selection.all().begin(), selection.all().end()) << std::endl;

  BOOST_CHECK_EQUAL(primclex.db<Supercell>().size(), 13);
  primclex.db<Supercell>().commit();
  DB::Selection<Supercell> selection2(primclex, "ALL");

  auto &qh2 = primclex.settings().query_handler<Supercell>();
  qh2.set_selected(selection2);
  auto &dict2 = qh2.dict();

  std::vector<std::string> all_col2 = {"name", "selected", "scel_size"};
  auto formatter2 = dict2.parse(all_col2);

  std::stringstream ss2;
  ss2 << formatter2(selection2.all().begin(), selection2.all().end()) << std::endl;

  BOOST_CHECK_EQUAL(primclex.db<PrimPeriodicDiffTransOrbit>().size(), 234);
  primclex.db<PrimPeriodicDiffTransOrbit>().commit();
  DB::Selection<PrimPeriodicDiffTransOrbit> selection3(primclex, "ALL");

  auto &qh3 = primclex.settings().query_handler<PrimPeriodicDiffTransOrbit>();
  qh3.set_selected(selection3);
  auto &dict3 = qh3.dict();

  std::vector<std::string> all_col3 = {"name", "selected"};
  auto formatter3 = dict3.parse(all_col3);

  std::stringstream ss3;
  ss3 << formatter3(selection3.all().begin(), selection3.all().end()) << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()




