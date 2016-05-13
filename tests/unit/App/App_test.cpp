#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
///   the command line executable

/// What is being used to test it:
#include <boost/filesystem.hpp>

#include "casm/CASM_classes.hh"
#include "Common.hh"

using namespace CASM;

struct Checks {

  Checks(std::string _type,
         std::string _no_proj_command,
         std::string _help_command) :
    type(_type),
    no_proj_command(_no_proj_command),
    help_command(_help_command) {}

  std::string type;
  std::string no_proj_command;
  std::string help_command;
};

BOOST_AUTO_TEST_SUITE(AppTest)

BOOST_AUTO_TEST_CASE(ProjectCommands) {

  Popen p;

  // checks for what happens when running casm from a non-project directory
  std::vector<Checks> command = {
    Checks("composition", "-d", "-h"),
    Checks("sym", "", "-h"),
    Checks("bset", "", "-h"),
    Checks("enum", "--supercells --max 1", "-h"),
    Checks("select", "", "-h"),
    Checks("query", "", "-h")
  };

  // check help doesn't need to be in a project
  p.popen("casm init -h");
  BOOST_CHECK_EQUAL(p.exit_code(), 0);

  for(auto it = command.begin(); it != command.end(); ++it) {
    p.popen("casm " + it->type + " " + it->no_proj_command);
    BOOST_CHECK_EQUAL(p.exit_code(), 3);

    p.popen("casm " + it->type + " " + it->help_command);
    BOOST_CHECK_EQUAL(p.exit_code(), 0);
  }

  // checks of 'casm X' commands for several projects
  std::vector<std::unique_ptr<test::Proj> > proj;
  proj.push_back(notstd::make_unique<test::FCCTernaryProj>());
  proj.push_back(notstd::make_unique<test::ZrOProj>());

  for(auto proj_it = proj.begin(); proj_it != proj.end(); ++proj_it) {
    (*proj_it)->check_init();
    (*proj_it)->check_symmetry();
    (*proj_it)->check_composition();
    (*proj_it)->check_bset();
    (*proj_it)->check_enum();
    (*proj_it)->check_select();
    (*proj_it)->check_query();
    rm_project(**proj_it);
  }

}

BOOST_AUTO_TEST_SUITE_END()
