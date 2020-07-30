#ifndef CASMtest_ZrOProj
#define CASMtest_ZrOProj

#include "gtest/gtest.h"
#include <boost/filesystem.hpp>

#include "crystallography/TestStructures.hh"
#include "Proj.hh"
#include "casm/casm_io/Log.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell_impl.hh"
#include "casm/clex/Configuration_impl.hh"
#include "casm/database/Database.hh"

using namespace CASM;

namespace test {

  class ZrOProj : public Proj {

  public:

    ZrOProj() :
      Proj(proj_dir(autotools::abs_srcdir() + "/tests/unit/test_projects/ZrO"),
           ZrO_prim(),
           "ZrO",
           "HCP Zr with octahedral interstitial O") {}

    /// \brief Check symmetry
    void check_symmetry() override {
      _check_symmetry(24, 12, 24, 12, 24, 12, "D6h", "D6h");
    }

    jsonParser bspecs() const {

      std::string str = R"({
"basis_function_specs" : {
  "dof_specs": {
    "occ": {
      "site_basis_functions" : "occupation"
    }
  }
},
"cluster_specs": {
  "method": "periodic_max_length",
  "params": {
    "orbit_branch_specs" : {
      "2" : {"max_length" : 9.0},
      "3" : {"max_length" : 7.0},
      "4" : {"max_length" : 6.0}
    },
    "orbit_specs" : [
      {
        "coordinate_mode" : "Integral",
        "prototype" : [
          [ 2, 0, 0, 0 ],
          [ 2, 3, 0, 0 ]
        ],
        "include_subclusters" : false
      }
    ]
  }
}
})";

      return jsonParser::parse(str);

    }

    void check_composition() override {
      std::vector<std::string> axes = {
        R"(\s+0\s+Zr\(2\)O\(2\)\s+Zr\(2\)Va\(2\)\s+Zr\(2\)Va\(2a\)O\(2-2a\))",
        R"(\s+1\s+Zr\(2\)Va\(2\)\s+Zr\(2\)O\(2\)\s+Zr\(2\)Va\(2-2a\)O\(2a\))"
      };

      _check_composition_axes(axes.begin(), axes.end());
    }

    //TODO: This code has been copied and pasted... there's probably a lot of it lurking around
    /// \brief Uses bspecs() and checks that 5 branches are generated,
    ///        and that --orbits, --clusters, and --functions run without error.
    void check_bset() override {

      // check for success with a valid bspecs
      bspecs().write(dir / "basis_sets" / "bset.default" / "bspecs.json");

      m_p.popen(cd_and() + autotools::abs_ccasm_path() + " bset -u");
      EXPECT_EQ(m_p.exit_code(), 0) << m_p.gets();

      EXPECT_EQ(boost::regex_search(m_p.gets(), m_match, boost::regex(R"(Write:.*clust\.json)")), true) << m_p.gets();
      EXPECT_EQ(boost::regex_search(m_p.gets(), m_match, boost::regex(R"(Write:.*basis\.json)")), true) << m_p.gets();
      EXPECT_EQ(boost::regex_search(m_p.gets(), m_match, boost::regex(R"(Write:.*)" + title + R"(_Clexulator\.cc)")), true) << m_p.gets();

      EXPECT_EQ(true, fs::exists(m_dirs.clust(m_set->default_clex().bset))) << m_p.gets();
      EXPECT_EQ(true, fs::exists(m_dirs.clexulator_src(m_set->project_name(), m_set->default_clex().bset))) << m_p.gets();

      std::string str;

      // check that 4 branches are created (null - 4-pt)
      // check that --orbits, --clusters, --functions all execute
      //   (derived Proj would have to check the actual results)
      std::string pattern = R"(\*\* Branch [0-9]+ \*\*)";
      boost::regex re(pattern);

      std::vector<std::string> checks = {
        autotools::abs_ccasm_path() + " bset --orbits",
        autotools::abs_ccasm_path() + " bset --clusters",
        autotools::abs_ccasm_path() + " bset --functions"
      };

      for(auto it = checks.begin(); it != checks.end(); ++it) {
        m_p.popen(cd_and() + *it);
        str = m_p.gets();

        auto begin = boost::sregex_iterator(str.begin(), str.end(), re);
        auto end = boost::sregex_iterator();
        auto count = std::distance(begin, end);

        EXPECT_EQ(count, 5) << m_p.gets();
      }

      // check that you can't overwrite without using -f
      m_p.popen(cd_and() + autotools::abs_ccasm_path() + " bset -u");
      EXPECT_EQ(m_p.exit_code(), 6);

      m_p.popen(cd_and() + autotools::abs_ccasm_path() + " bset -uf");
      EXPECT_EQ(m_p.exit_code(), 0);

    }

    /// \brief Check "ccasm enum"
    void check_enum() override {

      {
        m_p.popen(cd_and() + autotools::abs_ccasm_path() + " enum --method ScelEnum --max 10");
        std::stringstream ss;
        Log log(ss);
        PrimClex primclex(dir, log);
        EXPECT_EQ(primclex.generic_db<Supercell>().size(), 147) << m_p.gets();
      }

      {
        m_p.popen(cd_and() + autotools::abs_ccasm_path() + " enum --method ConfigEnumAllOccupations --max 6");
        std::stringstream ss;
        Log log(ss);
        PrimClex primclex(dir, log);
        EXPECT_EQ(primclex.generic_db<Configuration>().size(), 5763) << m_p.gets();
      }
    }

  };

}

#endif
