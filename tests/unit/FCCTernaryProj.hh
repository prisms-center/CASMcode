#ifndef CASMtest_FCCTernaryProj
#define CASMtest_FCCTernaryProj

using namespace CASM;

namespace test {
  
  inline BasicStructure<Site> FCC_ternary_prim() {
    
    // lattice vectors as rows
    Eigen::Matrix3d lat;
    lat << 2.0, 2.0, 0.0,
           0.0, 2.0, 2.0,
           2.0, 0.0, 2.0;
    
    BasicStructure<Site> struc(Lattice(lat.transpose()));
    struc.title = "FCC_ternary";
    
    Molecule A = make_atom("A", struc.lattice());
    Molecule B = make_atom("B", struc.lattice());
    Molecule C = make_atom("C", struc.lattice());
    
    struc.basis.push_back( Site(Coordinate(Eigen::Vector3d::Zero(), struc.lattice()), {A, B, C}) );
    
    return struc;
    
  }

  class FCCTernaryProj : public Proj {
  
  public:
  
    FCCTernaryProj() :
      Proj(fs::absolute(fs::path("tests/unit/App/FCC_ternary")), 
           FCC_ternary_prim(), 
           "FCC_ternary", 
           "FCC Ternary with A, B, C occupation") {}
    
    jsonParser bspecs() const {

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
      m_p.popen(cd_and() + "casm bset -u");
      BOOST_CHECK_EQUAL_MESSAGE(m_p.exit_code(), 4, m_p.gets());
      
      // check for success with a valid bspecs
      bspecs().write(dir / "basis_sets" / "bset.default" / "bspecs.json");
      
      m_p.popen(cd_and() + "casm bset -u");
      BOOST_CHECK_EQUAL_MESSAGE(m_p.exit_code(), 0, m_p.gets());
      
      BOOST_CHECK_EQUAL(std::regex_search(m_p.gets(), m_match, std::regex(R"(Wrote.*eci\.in)")), true);
      BOOST_CHECK_EQUAL(std::regex_search(m_p.gets(), m_match, std::regex(R"(Wrote.*clust\.json)")), true);
      BOOST_CHECK_EQUAL(std::regex_search(m_p.gets(), m_match, std::regex(R"(Wrote.*)" + title + R"(_Clexulator\.cc)")), true);
      
      BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.eci_in(m_set.bset())));
      BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.clust(m_set.bset())));
      BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.clexulator_src(m_set.name(), m_set.bset())));
      
      std::string str;
      
      // check that 5 branches are created (null - 4-pt)
      // check that --orbits, --clusters, --functions all execute 
      //   (derived Proj would have to check the actual results)
      std::string pattern = R"(\*\* Branch [0-9]+ \*\*)";
      std::regex re(pattern);
      
      std::vector<std::string> checks = {
        "casm bset --orbits", 
        "casm bset --clusters", 
        "casm bset --functions"
      };
      
      for(auto it=checks.begin(); it!=checks.end(); ++it) {
        m_p.popen(cd_and() + *it);
        str = m_p.gets();
      
        auto begin = std::sregex_iterator(str.begin(), str.end(), re);
        auto end = std::sregex_iterator();
        auto count = std::distance(begin, end);
        
        BOOST_CHECK_EQUAL(count, 5);
      }
      
      // check that you can't overwrite without using -f
      m_p.popen(cd_and() + "casm bset -u");
      BOOST_CHECK_EQUAL(m_p.exit_code(), 6);
      
      m_p.popen(cd_and() + "casm bset -uf");
      BOOST_CHECK_EQUAL(m_p.exit_code(), 0);
      
    }
    
    /// \brief Check "casm enum"
    void check_enum() override {
      
      {
        m_p.popen(cd_and() + "casm enum --supercells --max 10");
        std::stringstream ss;
        PrimClex primclex(dir, ss);
        BOOST_CHECK_EQUAL_MESSAGE(primclex.get_supercell_list().size(), 87, m_p.gets());
      }
      
      {
        m_p.popen(cd_and() + "casm enum --configs --max 6");
        std::stringstream ss;
        PrimClex primclex(dir, ss);
        BOOST_CHECK_EQUAL_MESSAGE(std::distance(primclex.config_begin(), primclex.config_end()), 1081, m_p.gets());
      }
    }
  
  };

}

#endif