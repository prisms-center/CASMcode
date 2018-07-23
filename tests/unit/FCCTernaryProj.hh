#ifndef CASMtest_FCCTernaryProj
#define CASMtest_FCCTernaryProj

#include "Proj.hh"
#include "casm/CASM_global_definitions.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigIterator.hh"


using namespace CASM;


namespace test {

  inline BasicStructure<Site> FCC_ternary_prim() {

    // lattice vectors as cols
    Eigen::Matrix3d lat;
    lat << 0.0, 2.0, 2.0,
        2.0, 0.0, 2.0,
        2.0, 2.0, 0.0;

    BasicStructure<Site> struc {Lattice{lat}};
    struc.title = "FCC_ternary";

    Molecule A = make_atom("A", struc.lattice());
    Molecule B = make_atom("B", struc.lattice());
    Molecule C = make_atom("C", struc.lattice());

    struc.basis.push_back(Site(Coordinate(Eigen::Vector3d::Zero(), struc.lattice(), CART), {A, B, C}));

    return struc;

  }

  class FCCTernaryProj : public Proj {

  public:

    FCCTernaryProj() :
      //Use PID to get unique naming. Otherwise different tests might obliterate your directory mid testing if you run in parallel
      Proj(proj_dir("tests/unit/test_projects/FCC_ternary"),
           FCC_ternary_prim(),
           "FCC_ternary",
           "FCC Ternary with A, B, C occupation") {}

    static jsonParser bspecs() {

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
        R"(\s+0\s+A\s+B\s+C\s+A\(1-a-b\)B\(a\)C\(b\))",
        R"(\s+1\s+B\s+A\s+C\s+A\(a\)B\(1-a-b\)C\(b\))",
        R"(\s+2\s+C\s+A\s+B\s+A\(a\)B\(b\)C\(1-a-b\))"
      };

      _check_composition_axes(axes.begin(), axes.end());
    }

    /// \brief Uses bspecs() and checks that 5 branches are generated,
    ///        and that --orbits, --clusters, and --functions run without error.
    void check_bset() override {

      // check for failure with bspecs with invalid JSON
      fs::ofstream file(dir / "basis_sets" / "bset.default" / "bspecs.json");
      file << invalid_bspecs() << "\n";
      file.close();
      m_p.popen(cd_and() + "ccasm bset -u");
      BOOST_CHECK_MESSAGE(m_p.exit_code() == 4, m_p.gets());

      // check for success with a valid bspecs
      bspecs().write(dir / "basis_sets" / "bset.default" / "bspecs.json");

      m_p.popen(cd_and() + "ccasm bset -u");
      BOOST_CHECK_MESSAGE(m_p.exit_code() == 0, m_p.gets());

      BOOST_CHECK_MESSAGE(boost::regex_search(m_p.gets(), m_match, boost::regex(R"(write:.*clust\.json)")) == true, m_p.gets());
      BOOST_CHECK_MESSAGE(boost::regex_search(m_p.gets(), m_match, boost::regex(R"(write:.*basis\.json)")) == true, m_p.gets());
      BOOST_CHECK_MESSAGE(boost::regex_search(m_p.gets(), m_match, boost::regex(R"(write:.*)" + title + R"(_Clexulator\.cc)")) == true, m_p.gets());

      BOOST_CHECK_MESSAGE(true == fs::exists(m_dirs.clust(m_set->default_clex().bset)), m_p.gets());
      BOOST_CHECK_MESSAGE(true == fs::exists(m_dirs.clexulator_src(m_set->name(), m_set->default_clex().bset)), m_p.gets());

      std::string str;

      // check that 5 branches are created (null - 4-pt)
      // check that --orbits, --clusters, --functions all execute
      //   (derived Proj would have to check the actual results)
      std::string pattern = R"(\*\* Branch [0-9]+ \*\*)";
      boost::regex re(pattern);

      std::vector<std::string> checks = {
        "ccasm bset --orbits",
        "ccasm bset --clusters",
        "ccasm bset --functions"
      };

      for(auto it = checks.begin(); it != checks.end(); ++it) {
        m_p.popen(cd_and() + *it);
        str = m_p.gets();

        auto begin = boost::sregex_iterator(str.begin(), str.end(), re);
        auto end = boost::sregex_iterator();
        auto count = std::distance(begin, end);

        BOOST_CHECK_MESSAGE(count == 5, m_p.gets());
      }

      // check that you can't overwrite without using -f
      m_p.popen(cd_and() + "ccasm bset -u");
      BOOST_CHECK_EQUAL(m_p.exit_code(), 6);

      m_p.popen(cd_and() + "ccasm bset -uf");
      BOOST_CHECK_EQUAL(m_p.exit_code(), 0);

    }

    /// \brief Check "ccasm enum"
    void check_enum() override {

      {
        m_p.popen(cd_and() + "ccasm enum --method ScelEnum --max 10");
        std::stringstream ss;
        Log log(ss);
        PrimClex primclex(dir, log);
        BOOST_CHECK_MESSAGE(primclex.get_supercell_list().size() == 87, m_p.gets());
      }

      {
        m_p.popen(cd_and() + "ccasm enum --method ConfigEnumAllOccupations --max 6");
        std::stringstream ss;
        Log log(ss);
        PrimClex primclex(dir, log);
        BOOST_CHECK_MESSAGE(std::distance(primclex.config_begin(), primclex.config_end()) == 1081, m_p.gets());
      }
    }

  };

}

#endif
