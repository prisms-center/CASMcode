#include "gtest/gtest.h"
#include "autotools.hh"

/// What is being tested:
#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/clex/ConfigEnumEquivalents.hh"
#include "casm/clex/ScelEnumEquivalents.hh"

/// What is being used to test it:
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/ProjectBuilder.hh"
#include "casm/database/Database.hh"
#include "Common.hh"

using namespace CASM;
using namespace test;

TEST(ConfigEnumEquivalentsTest, Test1) {

  // tests ConfigEnumAllOccupations and ConfigEnumEquivalents

  // read test file
  fs::path test_cases_path(autotools::abs_srcdir() + "/tests/unit/clex/ConfigEnumEquivalents_test_cases.json");
  jsonParser tests(test_cases_path);
  fs::path test_proj_dir(autotools::abs_srcdir() + "/tests/unit/clex/test_proj");

  for(auto test_it = tests.begin(); test_it != tests.end(); ++test_it) {

    // input and expected output data
    jsonParser &j = *test_it;

    // if false: print calculated results if no test data; if true: suppress
    bool quiet = false;
    j.get_else(quiet, "quiet", false);

    EXPECT_TRUE(j.contains("title")) << "test case 'title' is required";
    EXPECT_TRUE(j.contains("prim")) << "test case 'prim' is required";
    EXPECT_TRUE(j.contains("max_vol")) << "test case 'max_vol' is required";

    // generate prim
    Structure prim(read_prim(j["prim"], TOL));

    // clean up test proj
    if(fs::exists(test_proj_dir / ".casm")) {
      fs::remove_all(test_proj_dir);
    }

    fs::create_directory(test_proj_dir);

    j["prim"].write(test_proj_dir / "prim.json");

    // build a project
    ProjectBuilder builder(test_proj_dir, j["title"].get<std::string>(), "formation_energy");
    builder.build();

    // read primclex
    PrimClex primclex(test_proj_dir, null_log());
    double tol = primclex.crystallography_tol();

    // generate supercells
    ScelEnumProps enum_props(1, j["max_vol"].get<int>() + 1);
    ScelEnumByProps scel_enum(primclex, enum_props);
    for(const auto &scel : scel_enum) {
    }

    // generate configurations
    std::map<Index, Index> prim_count;
    std::map<Index, Index> total_count;

    Index i = 0;
    for(const auto &scel : primclex.generic_db<Supercell>()) {

      // for each supercell, enumerate unique configurations

      //std::cout << "scel: " << scel.name() << std::endl;
      //std::cout << "max_occupation: " << scel.max_allowed_occupation() << std::endl;

      ConfigEnumAllOccupations enumerator(scel);
      total_count[i] = 0;
      prim_count[i] = 0;

      // add prim configurations from smaller supercells that tile scel
      Index j = 0;
      auto scel_j = primclex.generic_db<Supercell>().begin();
      for(; scel_j->volume() < scel.volume(); ++scel_j) {
        jsonParser json;
        ScelEnumEquivalents e(*scel_j);
        for(const auto &tunit : e) {
          if(is_supercell(scel.lattice(), tunit.lattice(), tol).first) {
            total_count[i] += prim_count[j];
          }
        }
        ++j;
      }

      for(auto &config : enumerator) {
        ConfigEnumEquivalents equiv_enum(config);
        for(auto &equiv : equiv_enum) {
          ++total_count[i];
          ++prim_count[i];
          (void) equiv;
        }
      }

      auto max_occ = scel.max_allowed_occupation();
      long prod = std::accumulate(max_occ.begin(), max_occ.end(), long {1}, [](long a, int b) {
        return a * (b + 1);
      });

      // check that the number of Configurations enumerated is equal to
      // the product of the number of allowed occupants on each site in the
      // Supercell
      EXPECT_EQ(total_count[i], prod);

      ++i;
    }

    // run checks:
    // ... add more here ...

    // clean up test proj
    if(fs::exists(test_proj_dir / ".casm")) {
      fs::remove_all(test_proj_dir);
    }
  }
}
