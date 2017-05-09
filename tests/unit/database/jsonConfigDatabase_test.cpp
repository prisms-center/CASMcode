#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/database/ConfigDatabase.hh"
#include "casm/database/json/jsonDatabase.hh"


/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/casm_io/stream_io/container.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(ConfigDatabase_Test)

BOOST_AUTO_TEST_CASE(Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  const Structure &prim(primclex.prim());
  primclex.settings().set_crystallography_tol(1e-5);
  BOOST_CHECK_EQUAL(fs::equivalent(proj.dir, primclex.dir().root_dir()), true);

  DB::jsonDatabase<Configuration> db_config(primclex);

  db_config.open();
  BOOST_CHECK_EQUAL(db_config.size(), 0);

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = prim.lattice().vectors();

  Supercell tscel(
    &primclex,
    Lattice(2.*a, 2.*b, c));
  const Supercell &scel = *tscel.insert().first;

  Configuration config(scel, jsonParser(), ConfigDoF({0, 0, 0, 0}));
  auto res = db_config.insert(config);
  BOOST_CHECK_EQUAL(db_config.size(), 1);
  BOOST_CHECK_EQUAL(db_config.begin()->id(), "0");

  db_config.erase(res.first);
  BOOST_CHECK_EQUAL(db_config.size(), 0);

  ConfigEnumAllOccupations enum_config(scel);
  for(const auto &config : enum_config) {
    db_config.insert(config);
  }
  db_config.commit();
  BOOST_CHECK_EQUAL(db_config.size(), 12);
  BOOST_CHECK_EQUAL(db_config.begin()->id(), "1");
  for(const auto &config : db_config) {
    //std::cout << "id: " << config.id() << "  occ: " << config.occupation() << std::endl;
    BOOST_CHECK_EQUAL(config.cache().contains("multiplicity"), false);
    BOOST_CHECK_EQUAL(config.multiplicity() != 0, true);
    BOOST_CHECK_EQUAL(config.cache().contains("multiplicity"), true);
  }
  db_config.commit();

  {
    auto next = db_config.begin();
    auto it = next++;
    auto end = db_config.end();
    for(; next != end; ++it, ++next) {
      BOOST_CHECK_EQUAL(*it < *next, true);
    }
  }

  DB::jsonDB::DirectoryStructure jsonDB_dir(primclex.dir().root_dir());
  //jsonParser file(jsonDB_dir.obj_list<Configuration>());
  //file.print(std::cout);
  //fs::ifstream file(jsonDB_dir.obj_list<Configuration>());
  //std::cout << file.rdbuf() << std::endl;

  db_config.close();
  BOOST_CHECK_EQUAL(db_config.size(), 0);

  db_config.open();
  BOOST_CHECK_EQUAL(db_config.size(), 12);
  BOOST_CHECK_EQUAL(db_config.begin()->id(), "1");
  for(const auto &config : db_config) {
    //std::cout << "id: " << config.id() << "  occ: " << config.occupation() << std::endl;
    BOOST_CHECK_EQUAL(config.cache().contains("multiplicity"), true);
  }
  {
    auto next = db_config.begin();
    auto it = next++;
    auto end = db_config.end();
    for(; next != end; ++it, ++next) {
      BOOST_CHECK_EQUAL(*it < *next, true);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
