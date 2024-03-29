#ifndef CASMtest_FCCTernaryProj
#define CASMtest_FCCTernaryProj

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "Common.hh"
#include "Proj.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/Configuration_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell_impl.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/Site.hh"
#include "casm/database/Database.hh"
#include "casm/global/definitions.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

using namespace CASM;

namespace test {

class FCCTernaryProj : public Proj {
 public:
  FCCTernaryProj()
      : Proj(FCC_ternary_prim(), "FCC_ternary",
             "FCC Ternary with A, B, C occupation") {}

  static jsonParser bspecs() {
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
      "2" : {"max_length" : 4.01},
      "3" : {"max_length" : 3.01}
    },
    "orbit_specs" : [
      {
        "coordinate_mode" : "Direct",
        "prototype" : [
          [ 0.000000000000, 0.000000000000, 0.000000000000 ],
          [ 1.000000000000, 0.000000000000, 0.000000000000 ],
          [ 2.000000000000, 0.000000000000, 0.000000000000 ],
          [ 3.000000000000, 0.000000000000, 0.000000000000 ]
        ],
        "include_subclusters" : true
      },
      {
        "coordinate_mode" : "Direct",
        "prototype" : [
          [ 0.000000000000, 0.000000000000, 0.000000000000 ],
          [ 0.000000000000, 1.000000000000, 0.000000000000 ],
          [ 0.000000000000, 0.000000000000, 1.000000000000 ],
          [ 1.000000000000, 1.000000000000, 1.000000000000 ]
        ],
        "include_subclusters" : true
      }
    ]
  }
}
})";

    return jsonParser::parse(str);
  }

  std::string invalid_bspecs() const {
    std::string str = R"({
"basis_functions" : {
"site_basis_functions" : "occupation"
},
"orbit_branch_specs" : {
"2" : {"max_length" : 4.01},
"3" : {"max_length" : 3.01}
},
"orbit_specs" : [
{
"coordinate_mode" : "Direct",
"prototype" : [
[ 0.000000000000, 0.000000000000, 0.000000000000 ],
[ 1.000000000000, 0.000000000000, 0.000000000000 ],
[ 2.000000000000, 0.000000000000, 0.000000000000 ],
[ 3.000000000000, 0.000000000000, 0.000000000000 ],
],
"include_subclusters" : true
},
{
"coordinate_mode" : "Direct",
"prototype" : [
[ 0.000000000000, 0.000000000000, 0.000000000000 ],
[ 0.000000000000, 1.000000000000, 0.000000000000 ],
[ 0.000000000000, 0.000000000000, 1.000000000000 ],
[ 1.000000000000, 1.000000000000, 1.000000000000 ]
],
"include_subclusters" : true
}
]
})";

    return str;
  }

  /// \brief Check symmetry
  void check_symmetry() override {
    _check_symmetry(48, 10, 48, 10, 48, 10, "Oh", "Oh");
  }
  void check_composition() override {
    std::vector<std::string> axes = {
        R"(\s+0\s+C\s+B\s+A\s+A\(b\)B\(a\)C\(1-a-b\))",
        R"(\s+1\s+B\s+C\s+A\s+A\(b\)B\(1-a-b\)C\(a\))",
        R"(\s+2\s+A\s+C\s+B\s+A\(1-a-b\)B\(b\)C\(a\))"};

    _check_composition_axes(axes.begin(), axes.end());
  }

  // TODO: This code has been copied and pasted... there's probably a lot of it
  // lurking around
  /// \brief Uses bspecs() and checks that 5 branches are generated,
  ///        and that --orbits, --clusters, and --functions run without error.
  void check_bset() override {
    // check for failure with bspecs with invalid JSON
    fs::ofstream file(dir / "basis_sets" / "bset.default" / "bspecs.json");
    file << invalid_bspecs() << "\n";
    file.close();
    m_p.popen(cd_and() + autotools::abs_ccasm_path() + " bset -u");
    EXPECT_EQ(m_p.exit_code(), 4) << m_p.gets();

    // check for success with a valid bspecs
    bspecs().write(dir / "basis_sets" / "bset.default" / "bspecs.json");

    m_p.popen(cd_and() + autotools::abs_ccasm_path() + " bset -u");
    EXPECT_EQ(m_p.exit_code(), 0) << m_p.gets();

    EXPECT_EQ(true, fs::exists(m_dirs.clust(m_set->default_clex().bset)))
        << m_p.gets();
    EXPECT_EQ(true, fs::exists(m_dirs.basis(m_set->default_clex().bset)))
        << m_p.gets();
    EXPECT_EQ(true, fs::exists(m_dirs.clexulator_src(
                        m_set->project_name(), m_set->default_clex().bset)))
        << m_p.gets();

    std::string str;

    // check that 5 branches are created (null - 4-pt)
    // check that --orbits, --clusters, --functions all execute
    //   (derived Proj would have to check the actual results)
    std::string pattern = R"(\*\* Branch [0-9]+ \*\*)";
    boost::regex re(pattern);

    std::vector<std::string> checks = {
        autotools::abs_ccasm_path() + " bset --orbits",
        autotools::abs_ccasm_path() + " bset --clusters",
        autotools::abs_ccasm_path() + " bset --functions"};

    for (auto it = checks.begin(); it != checks.end(); ++it) {
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
      m_p.popen(cd_and() + autotools::abs_ccasm_path() +
                " enum --method ScelEnum --max 10");
      // ScopedStringStreamLogging logging;
      PrimClex primclex(dir);
      EXPECT_EQ(primclex.generic_db<Supercell>().size(), 87) << m_p.gets();
    }

    {
      m_p.popen(cd_and() + autotools::abs_ccasm_path() +
                " enum --method ConfigEnumAllOccupations --max 6");
      ScopedStringStreamLogging logging;
      PrimClex primclex(dir);
      EXPECT_EQ(primclex.generic_db<Configuration>().size(), 1081)
          << m_p.gets();
    }
  }
};

}  // namespace test

#endif
