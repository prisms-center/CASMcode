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

using namespace CASM;

BOOST_AUTO_TEST_SUITE(DataFormatterDictionaryTest)

BOOST_AUTO_TEST_CASE(Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();
  proj.check_composition();

  PrimClex primclex(proj.dir, null_log());
  const Structure &prim(primclex.prim());
  primclex.settings().set_crystallography_tol(1e-5);

  ScelEnumByProps enum_scel(primclex, ScelEnumProps(1, 5));
  ConfigEnumAllOccupations::run(primclex, enum_scel.begin(), enum_scel.end());
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
}

BOOST_AUTO_TEST_SUITE_END()




