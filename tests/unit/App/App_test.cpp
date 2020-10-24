#include "gtest/gtest.h"
#include "autotools.hh"

/// What is being tested:
///   the command line executable

/// What is being used to test it:
#include <boost/filesystem.hpp>

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "ZrOProj.hh"
#include "casm/system/Popen.hh"
#include "casm/casm_io/Log.hh"

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


TEST(AppTest, ProjectCommands) {

  ScopedNullLogging logging;
  Popen p;
  Log &log = CASM::log();

  // checks for what happens when running casm from a non-project directory
  std::vector<Checks> command = {
    Checks("composition", "-d", "-h"),
    Checks("sym", "", "-h"),
    Checks("bset", "", "-h"),
    Checks("enum", "", "-h"),
    Checks("select", "", "-h"),
    Checks("query", "", "-h")
  };

  // check help doesn't need to be in a project
  p.popen(autotools::abs_ccasm_path() + " init -h");
  EXPECT_EQ(p.exit_code(), 0) << p.gets();

  for(auto it = command.begin(); it != command.end(); ++it) {
    p.popen(autotools::abs_ccasm_path() + " " + it->type + " " + it->no_proj_command);
    EXPECT_EQ(p.exit_code(), 3) << p.gets();

    p.popen(autotools::abs_ccasm_path() + " " + it->type + " " + it->help_command);
    EXPECT_EQ(p.exit_code(), 0) << p.gets();
  }

  // checks of 'ccasm X' commands for several projects
  std::vector<std::unique_ptr<test::Proj> > proj;
  proj.push_back(notstd::make_unique<test::FCCTernaryProj>());
  proj.push_back(notstd::make_unique<test::ZrOProj>());

  for(auto proj_it = proj.begin(); proj_it != proj.end(); ++proj_it) {
    log.custom<Log::standard>("Test project: " + (*proj_it)->title);
    log << "root: " << (*proj_it)->dir << std::endl;

    log << "testing 'ccasm init'" << std::endl;
    (*proj_it)->check_init();

    log << "testing 'ccasm sym'" << std::endl;
    (*proj_it)->check_symmetry();

    log << "testing 'ccasm composition'" << std::endl;
    (*proj_it)->check_composition();

    log << "testing 'ccasm bset'" << std::endl;
    (*proj_it)->check_bset();

    log << "testing 'ccasm enum'" << std::endl;
    (*proj_it)->check_enum();

    log << "testing 'ccasm select'" << std::endl;
    (*proj_it)->check_select();

    log << "testing 'ccasm query'" << std::endl;
    (*proj_it)->check_query();

    log << "done" << std::endl;
    log << std::endl;

  }

  log << "delete test projects" << std::endl;

  log << "delete proj[0]:" << std::endl;
  proj[0].reset();
  log << "  done" << std::endl;

  log << "delete proj[1]:" << std::endl;
  proj[1].reset();
  log << "  done" << std::endl;

  log << "leaving test cast ProjectCommands" << std::endl;
}
