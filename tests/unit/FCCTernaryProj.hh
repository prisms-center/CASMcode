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