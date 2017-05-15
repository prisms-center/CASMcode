#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/database/json/jsonPropertiesDatabase.hh"


/// What is being used to test it:

#include <boost/filesystem.hpp>
#include "Common.hh"
#include "ZrOProj.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(jsonPropertiesDatabase_Test)

BOOST_AUTO_TEST_CASE(Test1) {

  test::ZrOProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  const Structure &prim(primclex.prim());
  primclex.settings().set_crystallography_tol(1e-5);
  BOOST_CHECK_EQUAL(fs::equivalent(proj.dir, primclex.dir().root_dir()), true);

  fs::path loc("tests/unit/database/config_props.json");
  DB::jsonPropertiesDatabase db_props(primclex, loc);
  BOOST_CHECK_EQUAL(1, 1);

  db_props.open();
  BOOST_CHECK_EQUAL(1, 1);
  BOOST_CHECK_EQUAL(db_props.empty(), true);
  BOOST_CHECK_EQUAL(db_props.size(), 0);

  DB::MappedProperties props;
  props.from = "from/0";
  props.to = "to/0";
  props.unmapped = jsonParser::parse(std::string(R"({"test":"data"})"));
  props.mapped = jsonParser::parse(std::string(R"({"relaxed_energy":0.1})"));

  auto res = db_props.insert(props);
  BOOST_CHECK_EQUAL(res.second, true);
  BOOST_CHECK_EQUAL(db_props.size(), 1);
  BOOST_CHECK_EQUAL(db_props.find_via_from("from/0")->to, "to/0");
  BOOST_CHECK_EQUAL(db_props.find_via_to("to/0")->from, "from/0");
  BOOST_CHECK_EQUAL(db_props.relaxed_from_all("to/0").size(), 1);
  BOOST_CHECK_EQUAL(db_props.relaxed_to("from/0"), "to/0");
  BOOST_CHECK_EQUAL(db_props.relaxed_from("to/0"), "from/0");
  BOOST_CHECK_EQUAL(db_props.best_score("to/0"), 0.1);
  BOOST_CHECK_EQUAL(db_props.score("from/0"), 0.1);

  props.from = "from/1";
  props.to = "to/1";
  props.unmapped = jsonParser::parse(std::string(R"({"test":"data"})"));
  props.mapped = jsonParser::parse(std::string(R"({"relaxed_energy":0.2})"));

  res = db_props.insert(props);
  BOOST_CHECK_EQUAL(res.second, true);
  BOOST_CHECK_EQUAL(db_props.size(), 2);
  BOOST_CHECK_EQUAL(db_props.find_via_from("from/1")->to, "to/1");
  BOOST_CHECK_EQUAL(db_props.find_via_to("to/1")->from, "from/1");
  BOOST_CHECK_EQUAL(db_props.relaxed_from_all("to/1").size(), 1);
  BOOST_CHECK_EQUAL(db_props.relaxed_to("from/1"), "to/1");
  BOOST_CHECK_EQUAL(db_props.relaxed_from("to/1"), "from/1");
  BOOST_CHECK_EQUAL(almost_equal(db_props.best_score("to/1"), 0.2), true);
  BOOST_CHECK_EQUAL(almost_equal(db_props.score("from/1"), 0.2), true);

  props.from = "from/2";
  props.to = "to/1";
  props.unmapped = jsonParser::parse(std::string(R"({"test":"data"})"));
  props.mapped = jsonParser::parse(std::string(R"({"relaxed_energy":0.3})"));

  res = db_props.insert(props);
  BOOST_CHECK_EQUAL(res.second, true);
  BOOST_CHECK_EQUAL(db_props.size(), 3);
  BOOST_CHECK_EQUAL(db_props.find_via_from("from/2")->to, "to/1");
  BOOST_CHECK_EQUAL(db_props.find_via_to("to/1")->from, "from/1");
  BOOST_CHECK_EQUAL(db_props.relaxed_from_all("to/1").size(), 2);
  BOOST_CHECK_EQUAL(db_props.relaxed_to("from/2"), "to/1");
  BOOST_CHECK_EQUAL(db_props.relaxed_from("to/1"), "from/1");
  BOOST_CHECK_EQUAL(almost_equal(db_props.best_score("to/1"), 0.2), true);
  BOOST_CHECK_EQUAL(almost_equal(db_props.score("from/2"), 0.3), true);

  auto it = db_props.find_via_to("to/1");
  db_props.erase(it);
  BOOST_CHECK_EQUAL(db_props.size(), 2);
  BOOST_CHECK_EQUAL(db_props.find_via_from("from/1") == db_props.end(), true);
  BOOST_CHECK_EQUAL(db_props.find_via_to("to/1")->from, "from/2");
  BOOST_CHECK_EQUAL(db_props.relaxed_from_all("to/1").size(), 1);
  BOOST_CHECK_EQUAL(db_props.relaxed_from("to/1"), "from/2");
  BOOST_CHECK_EQUAL(almost_equal(db_props.best_score("to/1"), 0.3), true);

  db_props.commit();
  BOOST_CHECK_EQUAL(1, 1);
  fs::ifstream file(loc);
  std::cout << file.rdbuf() << std::endl;

  db_props.close();
  BOOST_CHECK_EQUAL(1, 1);

  db_props.open();
  BOOST_CHECK_EQUAL(db_props.empty(), false);
  BOOST_CHECK_EQUAL(db_props.size(), 2);
  BOOST_CHECK_EQUAL(std::distance(db_props.begin(), db_props.end()), 2);

  db_props.close();
  BOOST_CHECK_EQUAL(1, 1);

  fs::remove(loc);
}

BOOST_AUTO_TEST_SUITE_END()
