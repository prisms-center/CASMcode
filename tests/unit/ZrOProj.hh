#ifndef CASMtest_ZrOProj
#define CASMtest_ZrOProj

#include "Proj.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigIterator.hh"

using namespace CASM;

namespace test {

  inline BasicStructure<Site> ZrO_prim() {

    // lattice vectors as rows
    Eigen::Matrix3d lat;
    lat << 3.233986860000, 0.000000000000, 0.000000000000,
        -1.616993430000, 2.800714770000, 0.000000000000,
        0.000000000000, 0.000000000000, 5.168678340000;

    BasicStructure<Site> struc(Lattice(lat.transpose()));
    struc.title = "ZrO";

    Molecule O = make_atom("O", struc.lattice());
    Molecule Zr = make_atom("Zr", struc.lattice());
    Molecule Va = make_vacancy(struc.lattice());

    struc.basis.push_back(Site(Coordinate(Eigen::Vector3d::Zero(), struc.lattice(), FRAC), {Zr}));
    struc.basis.push_back(Site(Coordinate(Eigen::Vector3d(2. / 3., 1. / 3., 1. / 2.), struc.lattice(), FRAC), {Zr}));
    struc.basis.push_back(Site(Coordinate(Eigen::Vector3d(1. / 3., 2. / 3., 1. / 4.), struc.lattice(), FRAC), {Va, O}));
    struc.basis.push_back(Site(Coordinate(Eigen::Vector3d(1. / 3., 2. / 3., 3. / 4.), struc.lattice(), FRAC), {Va, O}));

    return struc;
  }

  class ZrOProj : public Proj {

  public:

    ZrOProj() :
      //Use PID to get unique naming. Otherwise different tests might obliterate your directory mid testing if you run in parallel
      Proj(proj_dir("tests/unit/test_projects/ZrO"),
           ZrO_prim(),
           "ZrO",
           "HCP Zr with octahedral interstitial O") {}

    /// \brief Check symmetry
    void check_symmetry() override {
      _check_symmetry(24, 12, 24, 12, 24, 12, "D6h", "D6h");
    }

    jsonParser bspecs() const {

      std::string str = R"({
"basis_functions" : {
"site_basis_functions" : "occupation"
},
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
})";

      return jsonParser::parse(str);

    }

    void check_composition() override {
      std::vector<std::string> axes = {
        R"(\s+0\s+Zr\(2\)Va\(2\)\s+Zr\(2\)O\(2\)\s+Zr\(2\)Va\(2-2a\)O\(2a\))",
        R"(\s+1\s+Zr\(2\)O\(2\)\s+Zr\(2\)Va\(2\)\s+Zr\(2\)Va\(2a\)O\(2-2a\))"
      };

      _check_composition_axes(axes.begin(), axes.end());
    }

    /// \brief Uses bspecs() and checks that 5 branches are generated,
    ///        and that --orbits, --clusters, and --functions run without error.
    void check_bset() override {

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

      // check that 4 branches are created (null - 4-pt)
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
        BOOST_CHECK_MESSAGE(primclex.get_supercell_list().size() == 147, m_p.gets());
      }

      {
        m_p.popen(cd_and() + "ccasm enum --method ConfigEnumAllOccupations --max 6");
        std::stringstream ss;
        Log log(ss);
        PrimClex primclex(dir, log);
        BOOST_CHECK_MESSAGE(std::distance(primclex.config_begin(), primclex.config_end()) == 5763, m_p.gets());
      }
    }

  };

}

#endif
