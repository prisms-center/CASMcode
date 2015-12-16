#ifndef CASM_unit_Common
#define CASM_unit_Common

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

// adapted from: http://stackoverflow.com/questions/26652904/boost-check-equal-and-dervatives-add-custom-message
//____________________________________________________________________________//

#define BOOST_TEST_REL_EQ_MESSAGE_EXTENSION(L, R, M, CMP, ICMP, CT)         \
    {                                                                       \
        auto _1(L);                                                         \
        auto _2(R);                                                         \
        std::stringstream ss;                                               \
        ss << "check " << BOOST_TEST_STRINGIZE(L) << " " << BOOST_TEST_STRINGIZE(CMP) \
           << " " << BOOST_TEST_STRINGIZE(R) << " failed [" << _1 << " " \
           << BOOST_TEST_STRINGIZE(ICMP) << " " << _2 << "] : \n---\n" << M \
           << "\n---\n";\
        BOOST_CHECK_IMPL( (_1 CMP _2), ss.str(), CT, CHECK_MSG );           \
    }                                                                       \
/**/

#define BOOST_CHECK_EQUAL_MESSAGE(L, R, M)      BOOST_TEST_REL_EQ_MESSAGE_EXTENSION(L, R, M, ==, !=, CHECK )
#define BOOST_WARN_EQUAL_MESSAGE(L, R, M)       BOOST_TEST_REL_EQ_MESSAGE_EXTENSION(L, R, M, ==, !=, WARN )
#define BOOST_REQUIRE_EQUAL_MESSAGE(L, R, M)    BOOST_TEST_REL_EQ_MESSAGE_EXTENSION(L, R, M, ==, !=, REQUIRE )

#include "casm/core"

using namespace CASM;

namespace test {
  
  class Proj {
  
  public:
  
    Proj(fs::path _proj_dir, 
             const BasicStructure<Site>& _prim, 
             std::string _title, 
             std::string _desc) :
      dir(_proj_dir),
      prim(_prim),
      title(_title),
      desc(_desc),
      m_dirs(dir) {}
    
    fs::path dir;
    BasicStructure<Site> prim;
    std::string title;
    std::string desc;
    
    std::string cd_and() const {
      return "cd " + dir.string() + "&& ";
    };
    
    virtual jsonParser bspecs() const;
    
    /// \brief Check project initialization
    virtual void check_init();
    
    /// \brief Check symmetry
    virtual void check_symmetry();
    
    /// \brief Check "casm composition"
    virtual void check_composition();
    
    /// \brief Check "casm bset"
    virtual void check_bset();
    
    /// \brief Check "casm enum"
    virtual void check_enum();
    
    /// \brief Check "casm select"
    virtual void check_select();
    
    /// \brief Check "casm query"
    virtual void check_query();
  
  
  protected:
    
    Popen m_p;
    std::smatch m_match;
    DirectoryStructure m_dirs;
    ProjectSettings m_set;
    
    void _check_symmetry(int lat_pg_op, int lat_pg_class,
                         int xtal_pg_op, int xtal_pg_class,
                         int fg_op, int fg_class,
                         std::string lat_pg_name, std::string xtal_pg_name);
    
    template<typename Iterator>
    void _check_composition_axes(Iterator begin, Iterator end);
    
  };

  /// \brief Build a CASM project at 'proj_dir' using the prim
  ///
  /// Typically, proj_dir.filename() == title 
  void make_project(const Proj& proj);
  
  /// \brief Remove a CASM project, checking first that there is a '.casm' dir
  ///
  /// Be careful! This does a recursive remove of the entire proj_dir!
  void rm_project(const Proj& proj);
  
  /// \brief Check some aspects of a SymGroup json, including the expected
  ///        number of conjugacy classes and operations
  void check_symgroup(const jsonParser& json, int N_op, int N_class);
  
}

#include "FCCTernaryProj.hh"

namespace test {
  
  template<typename Iterator>
  void Proj::_check_composition_axes(Iterator begin, Iterator end) {
    
    m_p.popen(cd_and() + "casm composition --select 0");
  
    for(auto it=begin; it!=end; ++it) {
      BOOST_CHECK_EQUAL(std::regex_search(m_p.gets(), m_match, std::regex(*it)), true);
    }
    
    BOOST_CHECK_EQUAL(std::regex_search(m_p.gets(), m_match, std::regex(R"(Currently selected composition axes: 0)")), true);
  }
}

#endif