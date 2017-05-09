#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/database/Selection.hh"
#include "casm/database/DatabaseDefs.hh"


/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/ConfigEnumAllOccupations_impl.hh"
#include "casm/clex/ScelEnum.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(Selection_Test)

BOOST_AUTO_TEST_CASE(Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  const Structure &prim(primclex.prim());
  primclex.settings().set_crystallography_tol(1e-5);

  ScelEnumByProps enum_scel(primclex, ScelEnumProps(1, 5));
  ConfigEnumAllOccupations::run(primclex, enum_scel.begin(), enum_scel.end());
  BOOST_CHECK_EQUAL(primclex.db<Configuration>().size(), 126);

  DB::Selection<Configuration> selection(primclex);
  BOOST_CHECK_EQUAL(selection.size(), 126);
  BOOST_CHECK_EQUAL(selection.selected_size(), 0);

  auto &dict = primclex.settings().query_handler<Configuration>().dict();
  selection.set(dict, "lt(scel_size,3)");
  BOOST_CHECK_EQUAL(selection.selected_size(), 9);

}

BOOST_AUTO_TEST_SUITE_END()
