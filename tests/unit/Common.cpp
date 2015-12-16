#include "Common.hh"

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/app/ProjectBuilder.hh"

namespace test {

  BasicStructure<Site> ZrO_prim() {
    
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
    
    struc.basis.push_back( Site(Coordinate(Eigen::Vector3d::Zero(), struc.lattice()), {Zr}) );
    struc.basis.push_back( Site(Coordinate(Eigen::Vector3d(2./3., 1./3., 1./2.), struc.lattice()), {Zr}) );
    struc.basis.push_back( Site(Coordinate(Eigen::Vector3d(1./3., 2./3., 1./4.), struc.lattice()), {Va, O}) );
    struc.basis.push_back( Site(Coordinate(Eigen::Vector3d(1./3., 2./3., 3./4.), struc.lattice()), {Va, O}) );
    
    return struc;
  }

  /// \brief Build a CASM project at 'proj_dir/title' using the prim
  void make_project(const Proj& proj) {
    
    jsonParser json;
    write_prim(proj.prim, json, FRAC);
    json["description"] = proj.desc;
    
    fs::create_directory(proj.dir);
    
    json.write(proj.dir / "prim.json");
    
    // build a project
    ProjectBuilder builder(proj.dir, proj.title, "formation_energy");
    builder.build();
    
  }
  
  /// \brief Remove a CASM project, checking first that there is a '.casm' dir
  ///
  /// Be careful! This does a recursive remove of the entire proj_dir!
  void rm_project(const Proj& proj) {
    
    if(fs::exists(proj.dir/ ".casm")) {
      fs::remove_all(proj.dir);
    }
    
  }
  
  /// \brief Check some aspects of a SymGroup json
  void check_symgroup(const jsonParser& json, int N_op, int N_class) {
    BOOST_CHECK_EQUAL(json["character_table"].size(), N_class);
    BOOST_CHECK_EQUAL(json["conjugacy_class"].size(), N_class);
    
    BOOST_CHECK_EQUAL(json["symop"].size(), N_op);
    BOOST_CHECK_EQUAL(json["symop"][0]["type"].get<std::string>(), "identity");
    
    BOOST_CHECK_EQUAL(json["inverse"].size(), N_op);
    BOOST_CHECK_EQUAL(json["multiplication_table"].size(), N_op);
    for(auto i=0; i<json["multiplication_table"].size(); ++i) {
      BOOST_CHECK_EQUAL(json["multiplication_table"][i].size(), N_op);
    }
  }
  
  jsonParser Proj::bspecs() const {

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
  
  std::string Proj::invalid_bspecs() const {

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
  
  /// \brief Check that 'casm init' runs without error and expected files are created
  void Proj::check_init() {
      
      test::rm_project(*this);
    
      test::make_project(*this);
      
      m_set = ProjectSettings(dir);
    
      BOOST_CHECK_EQUAL(true, fs::exists(dir));
      
      // prim and settings
      BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.prim()));
      BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.project_settings()));
      
      // symmetry
      BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.crystal_point_group()));
      BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.factor_group()));
      BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.lattice_point_group()));
        
      // composition axes
      BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.composition_axes(m_set.calctype(), m_set.ref())));
    }
  
  /// \brief Check that 'casm sym' runs without error
  void Proj::check_symmetry() {
    m_p.popen(cd_and() + "casm sym");
    BOOST_CHECK_EQUAL(m_p.exit_code(), 0);
  }
  
  /// \brief Check number of symmetry operations and classes found
  void Proj::_check_symmetry(int lat_pg_op, int lat_pg_class,
                       int xtal_pg_op, int xtal_pg_class,
                       int fg_op, int fg_class,
                       std::string lat_pg_name, std::string xtal_pg_name) {
    
    
    check_symgroup(jsonParser(m_dirs.lattice_point_group()), lat_pg_op, lat_pg_class);
    check_symgroup(jsonParser(m_dirs.crystal_point_group()), xtal_pg_op, xtal_pg_class);
    check_symgroup(jsonParser(m_dirs.factor_group()), fg_op, fg_class);
    
    m_p.popen(cd_and() + "casm sym");
    
    BOOST_CHECK_EQUAL(std::regex_search(m_p.gets(), m_match, std::regex(R"(Lattice point group is:\s+)" + lat_pg_name)), true);
    BOOST_CHECK_EQUAL(std::regex_search(m_p.gets(), m_match, std::regex(R"(Crystal point group is:\s+)" + xtal_pg_name)), true);
  }
  
  /// \brief Default checks '-d' runs without error
  void Proj::check_composition() {
    m_p.popen(cd_and() + "casm composition -d");
    BOOST_CHECK_EQUAL(m_p.exit_code(), 0);
  }
    
  /// \brief Default uses Proj::bspecs and checks that 5 branches are generated,
  ///        and that --orbits, --clusters, and --functions run without error.
  void Proj::check_bset() {
    
    // check for failure with bspecs with invalid JSON
    fs::ofstream file(dir / "basis_sets" / "bset.default" / "bspecs.json");
    file << Proj::invalid_bspecs() << "\n";
    file.close();
    m_p.popen(cd_and() + "casm bset -u");
    BOOST_CHECK_EQUAL_MESSAGE(m_p.exit_code(), 4, m_p.gets());
    
    // check for success with a valid bspecs
    Proj::bspecs().write(dir / "basis_sets" / "bset.default" / "bspecs.json");
    
    m_p.popen(cd_and() + "casm bset -u");
    BOOST_CHECK_EQUAL(m_p.exit_code(), 0);
    
    BOOST_CHECK_EQUAL(std::regex_search(m_p.gets(), m_match, std::regex(R"(Wrote.*eci\.in)")), true);
    BOOST_CHECK_EQUAL(std::regex_search(m_p.gets(), m_match, std::regex(R"(Wrote.*clust\.json)")), true);
    BOOST_CHECK_EQUAL(std::regex_search(m_p.gets(), m_match, std::regex(R"(Wrote.*prim_nlist\.json)")), true);
    BOOST_CHECK_EQUAL(std::regex_search(m_p.gets(), m_match, std::regex(R"(Wrote.*)" + title + R"(_Clexulator\.cc)")), true);
    
    BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.eci_in(m_set.bset())));
    BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.clust(m_set.bset())));
    BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.prim_nlist(m_set.bset())));
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
  
  /// \brief Default checks that enumerating supercells and configurations can 
  /// be run for '--max 2' without error, but doesn't check results
  void Proj::check_enum() {
    m_p.popen(cd_and() + "casm enum --supercells --max 2");
    BOOST_CHECK_EQUAL(m_p.exit_code(), 0);
    m_p.popen(cd_and() + "casm enum --configs --max 2");
    BOOST_CHECK_EQUAL(m_p.exit_code(), 0);
  }
  
  /// \brief Default checks that several options run without error
  void Proj::check_select() {
    m_p.popen(cd_and() + "casm select --set-on");
    BOOST_CHECK_EQUAL(m_p.exit_code(), 0);
    
    m_p.popen(cd_and() + "casm select --set-off");
    BOOST_CHECK_EQUAL(m_p.exit_code(), 0);
  }
  
  /// \brief Default checks that several options run without error
  void Proj::check_query() {
    m_p.popen(cd_and() + "casm query --columns comp");
    BOOST_CHECK_EQUAL_MESSAGE(m_p.exit_code(), 0, m_p.gets());
    
  }
  
}

