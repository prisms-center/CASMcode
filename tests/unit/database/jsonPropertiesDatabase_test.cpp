#include "gtest/gtest.h"

/// What is being tested:
#include "casm/database/json/jsonPropertiesDatabase.hh"


/// What is being used to test it:

#include <boost/filesystem.hpp>
#include "Common.hh"
#include "ZrOProj.hh"

using namespace CASM;

TEST(jsonPropertiesDatabase_Test, Test1) {

  test::ZrOProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  //const Structure &prim(primclex.prim());
  primclex.settings().set_crystallography_tol(1e-5);
  EXPECT_EQ(fs::equivalent(proj.dir, primclex.dir().root_dir()), true);

  std::string calc_type("test");
  fs::path loc("tests/unit/database/config_props.json");
  DB::jsonPropertiesDatabase db_props(primclex, calc_type, loc);
  EXPECT_EQ(1, 1);

  db_props.open();
  EXPECT_EQ(1, 1);
  EXPECT_EQ(db_props.empty(), true);
  EXPECT_EQ(db_props.size(), 0);

  MappedProperties props;
  props.from = "from/0";
  props.to = "to/0";
  props.scalar("relaxed_energy") = 0.1;
  props.site["test"] = Eigen::MatrixXd::Ones(3, 3);

  auto res = db_props.insert(props);
  EXPECT_EQ(res.second, true);
  EXPECT_EQ(db_props.size(), 1);
  EXPECT_EQ(db_props.find_via_from("from/0")->to, "to/0");
  EXPECT_EQ(db_props.find_via_to("to/0")->from, "from/0");
  EXPECT_EQ(db_props.relaxed_from_all("to/0").size(), 1);
  EXPECT_EQ(db_props.relaxed_to("from/0"), "to/0");
  EXPECT_EQ(db_props.relaxed_from("to/0"), "from/0");
  EXPECT_EQ(db_props.best_score("to/0"), 0.1);
  EXPECT_EQ(db_props.score("from/0"), 0.1);

  props.from = "from/1";
  props.to = "to/1";
  props.scalar("relaxed_energy") = 0.2;
  props.site["test"] = Eigen::MatrixXd::Ones(3, 3);

  res = db_props.insert(props);
  EXPECT_EQ(res.second, true);
  EXPECT_EQ(db_props.size(), 2);
  EXPECT_EQ(db_props.find_via_from("from/1")->to, "to/1");
  EXPECT_EQ(db_props.find_via_to("to/1")->from, "from/1");
  EXPECT_EQ(db_props.relaxed_from_all("to/1").size(), 1);
  EXPECT_EQ(db_props.relaxed_to("from/1"), "to/1");
  EXPECT_EQ(db_props.relaxed_from("to/1"), "from/1");
  EXPECT_EQ(almost_equal(db_props.best_score("to/1"), 0.2), true);
  EXPECT_EQ(almost_equal(db_props.score("from/1"), 0.2), true);

  props.from = "from/2";
  props.to = "to/1";
  props.scalar("relaxed_energy") = 0.3;
  props.site["test"] = Eigen::MatrixXd::Ones(3, 3);

  res = db_props.insert(props);
  EXPECT_EQ(res.second, true);
  EXPECT_EQ(db_props.size(), 3);
  EXPECT_EQ(db_props.find_via_from("from/2")->to, "to/1");
  EXPECT_EQ(db_props.find_via_to("to/1")->from, "from/1");
  EXPECT_EQ(db_props.relaxed_from_all("to/1").size(), 2);
  EXPECT_EQ(db_props.relaxed_to("from/2"), "to/1");
  EXPECT_EQ(db_props.relaxed_from("to/1"), "from/1");
  EXPECT_EQ(almost_equal(db_props.best_score("to/1"), 0.2), true);
  EXPECT_EQ(almost_equal(db_props.score("from/2"), 0.3), true);

  auto it = db_props.find_via_to("to/1");
  db_props.erase(it);
  EXPECT_EQ(db_props.size(), 2);
  EXPECT_EQ(db_props.find_via_from("from/1") == db_props.end(), true);
  EXPECT_EQ(db_props.find_via_to("to/1")->from, "from/2");
  EXPECT_EQ(db_props.relaxed_from_all("to/1").size(), 1);
  EXPECT_EQ(db_props.relaxed_from("to/1"), "from/2");
  EXPECT_EQ(almost_equal(db_props.best_score("to/1"), 0.3), true);

  db_props.commit();
  EXPECT_EQ(1, 1);
  //fs::ifstream file(loc);
  //std::cout << file.rdbuf() << std::endl;

  db_props.close();
  EXPECT_EQ(1, 1);

  db_props.open();
  EXPECT_EQ(db_props.empty(), false);
  EXPECT_EQ(db_props.size(), 2);
  EXPECT_EQ(std::distance(db_props.begin(), db_props.end()), 2);

  db_props.close();
  EXPECT_EQ(1, 1);

  fs::remove(loc);
}
