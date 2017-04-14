#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/ScelEnum.hh"
//#include "casm/clex/ScelEnum_impl.hh"

/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/app/casm_functions.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(ScelEnumTest)

BOOST_AUTO_TEST_CASE(Test1) {

  test::ZrOProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  primclex.settings().set_crystallography_tol(1e-5);

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.get_prim().lattice().vectors();

  std::vector<std::string> m_names;

  // -- Test ScelEnumByProps --------------------
  {
    ScelEnumProps enum_props(1, 10);
    ScelEnumByProps e(primclex, enum_props);

    BOOST_CHECK_EQUAL(e.name(), "ScelEnumByProps");

    auto it = e.begin();
    BOOST_CHECK(true);
    BOOST_CHECK_EQUAL(it.name(), "ScelEnumByProps");

    auto end = e.end();
    BOOST_CHECK(true);

    Index count = 0;
    for(; it != end; ++it, ++count) {
      m_names.push_back(it->get_name());

      Lattice canon_check = canonical_equivalent_lattice(
                              it->get_real_super_lattice(),
                              primclex.get_prim().point_group(),
                              primclex.crystallography_tol());

      bool check = almost_equal(
                     it->get_real_super_lattice().lat_column_mat(),
                     canon_check.lat_column_mat(),
                     primclex.crystallography_tol());

      if(!check) {
        std::cout << "superlat: \n" << it->get_real_super_lattice().lat_column_mat() << std::endl;
        std::cout << "canon_check: \n" << canon_check.lat_column_mat() << std::endl;
      }

      BOOST_CHECK_EQUAL(check, true);

      BOOST_CHECK_EQUAL(it->is_canonical(), true);

    }
    BOOST_CHECK_EQUAL(count, 114);
    BOOST_CHECK(it == end);
  }

  // -- use results to Test ScelEnumByName --------------------
  {
    ScelEnumByNameT<false> e(primclex, m_names.begin(), m_names.end());
    BOOST_CHECK_EQUAL(e.name(), "ScelEnumByName");

    auto it = e.begin();
    BOOST_CHECK(true);
    BOOST_CHECK_EQUAL(it.name(), "ScelEnumByName");

    auto end = e.end();
    BOOST_CHECK(true);

    Index count = 0;
    for(; it != end; ++it, ++count) {
      //std::cout << it->get_name() << std::endl;
    }
    BOOST_CHECK_EQUAL(count, 114);
    BOOST_CHECK(it == end);
  }

}

BOOST_AUTO_TEST_CASE(Test2) {

  // create a project
  test::FCCTernaryProj proj;
  proj.check_init();

  // in case you want to see what's happening
  OStringStreamLog ss_log;
  OStringStreamLog ss_debug_log;
  OStringStreamLog ss_err_log;

  // construct PrimClex
  PrimClex primclex(proj.dir, Logging(ss_log, ss_debug_log, ss_err_log));

  auto exec = [&](const std::string & args) {
    CommandArgs cmdargs(args, &primclex, proj.dir, ss_log, ss_err_log);
    int code = casm_api(cmdargs);
    //std::cout << ss_log.ss().str() << std::endl;
    //std::cout << ss_err_log.ss().str() << std::endl;
    return code;
  };

  BOOST_CHECK_EQUAL(exec("casm enum -h"), 0);
  BOOST_CHECK_EQUAL(exec("casm enum --method ScelEnum --max 4"), 0);
  BOOST_CHECK_EQUAL(exec("casm enum --method ConfigEnumAllOccupations --all"), 0);
  BOOST_CHECK_EQUAL(exec("casm enum --method ScelEnum --max 8"), 0);
  BOOST_CHECK_EQUAL(exec("casm enum --method ConfigEnumAllOccupations --max 6 -i '{\"existing_only\":true}'"), 0);

}

BOOST_AUTO_TEST_SUITE_END()
