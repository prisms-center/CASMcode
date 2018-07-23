#ifndef CASM_UNIT_PROJ
#define CASM_UNIT_PROJ

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "casm/CASM_global_definitions.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/system/Popen.hh"

using namespace CASM;

namespace test {

  class Proj {

  public:

    Proj(fs::path _proj_dir,
         const BasicStructure<Site> &_prim,
         std::string _title,
         std::string _desc) :
      dir(_proj_dir),
      prim(_prim),
      title(_title),
      desc(_desc),
      m_dirs(dir) {}

    virtual ~Proj() {}

    fs::path dir;
    BasicStructure<Site> prim;
    std::string title;
    std::string desc;

    std::string cd_and() const {
      return "cd " + dir.string() + "&& ";
    };

    /// \brief Check project initialization
    virtual void check_init();

    /// \brief Check symmetry
    virtual void check_symmetry();

    /// \brief Check "casm composition"
    virtual void check_composition();

    /// \brief Check "casm bset"
    virtual void check_bset() = 0;

    /// \brief Check "casm enum"
    virtual void check_enum();

    /// \brief Check "casm select"
    virtual void check_select();

    /// \brief Check "casm query"
    virtual void check_query();

    /// \brief Build a CASM project at 'proj_dir' using the prim
    ///
    /// Typically, proj_dir.filename() == title
    void make();

    static std::vector<fs::path> &directory_list;

  protected:

    Popen m_p;
    boost::smatch m_match;
    DirectoryStructure m_dirs;
    notstd::cloneable_ptr<ProjectSettings> m_set;

    void _check_symmetry(int lat_pg_op, int lat_pg_class,
                         int xtal_pg_op, int xtal_pg_class,
                         int fg_op, int fg_class,
                         std::string lat_pg_name, std::string xtal_pg_name);

    template<typename Iterator>
    void _check_composition_axes(Iterator begin, Iterator end);

  };

  /// \brief Create a new project directory, appending ".(#)" to ensure
  /// it is a new project
  fs::path proj_dir(fs::path init);

  /// \brief Check some aspects of a SymGroup json, including the expected
  ///        number of conjugacy classes and operations
  void check_symgroup(const jsonParser &json, int N_op, int N_class);

}

namespace test {

  template<typename Iterator>
  void Proj::_check_composition_axes(Iterator begin, Iterator end) {

    m_p.popen(cd_and() + "ccasm composition --select 0");

    for(auto it = begin; it != end; ++it) {
      BOOST_CHECK_MESSAGE(boost::regex_search(m_p.gets(), m_match, boost::regex(*it)) == true, m_p.gets());
    }

    BOOST_CHECK_MESSAGE(
      boost::regex_search(m_p.gets(), m_match, boost::regex(R"(Currently selected composition axes: 0)")) ==
      true,
      m_p.gets()
    );
  }
}

#endif
