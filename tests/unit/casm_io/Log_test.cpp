#include "gtest/gtest.h"

#include "casm/casm_io/Log.hh"

using namespace CASM;

TEST(LogTest, BasicUsage) {

  // The default logs are initially available from CASM::log() and CASM::err_log():

  log() << "standard out" << std::endl;
  err_log() << "standard err" << std::endl;

}

TEST(LogTest, Sections) {

  // The log can be separated into sections and given a verbosity level to control which sections
  // are printed. The verbosity level is an integer. For convenience, some levels are named:
  // - Log::none == 0
  // - Log::quiet == 5
  // - Log::standard == 10
  // - Log::verbose == 20
  // - Log::debug == 100
  //
  // Any section that has a equal or lesser verbosity than the current verbosity_level will be
  // printed. Any section that has a higher verbosity level will be skipped.

  EXPECT_EQ(log().verbosity(), Log::standard);
  EXPECT_EQ(log().verbosity(), 10);

  log().begin_section<Log::standard>();
  log() << "standard verbosity level section (printed by default)" << std::endl;
  log().end_section();

  log().begin_section<Log::verbose>();
  log() << "verbose section (not printed by default)" << std::endl;
  log().end_section();

  // set verbosity level
  log().set_verbosity(Log::verbose);

  log().begin_section<Log::verbose>();
  log() << "verbose section (now it is printed)" << std::endl;
  log().end_section();

  // Sections may be nested. Ending a section ends the last section that was begun, reverting
  // the verbosity level so that the previous section's verbosity level becomes controlling again.
  // This features allows an individual function to control its own verbosity level.
  //
  log().set_verbosity(Log::standard);
  log().begin_section<Log::standard>(); // begin standard section
  log() << "standard section (is printed)" << std::endl;
  {
    log().begin_section<Log::verbose>(); // begin verbose section
    log() << "verbose section (is not printed)" << std::endl;
    log().end_section(); // end verbose section
  }
  log() << "standard section again (is printed)" << std::endl;
  log().end_section(); // end standard section

}

TEST(LogTest, NamedSections) {

  // Named sections:
  // - automatically end the last section and begin a new section
  // - print a section header with a standardized format

  log().begin<Log::standard>("Standard verbosity section");
  log() << "standard verbosity level section (printed by default)" << std::endl;

  log().begin<Log::verbose>("Verbose section");
  log() << "verbose section (not printed by default)" << std::endl;

  // set verbosity level
  log().set_verbosity(Log::verbose);

  log().begin<Log::verbose>("Verbose section");
  log() << "verbose section (now it is printed)" << std::endl;

  // To nest a named section, use the "subsection" and "end_section" methods.
  log().set_verbosity(Log::standard);

  log().begin<Log::standard>("Standard section"); // begin standard section
  log() << "standard section (is printed)" << std::endl;
  {
    log().subsection().begin<Log::verbose>("Verbose subsection"); // begin verbose section
    log() << "verbose section (is not printed)" << std::endl;
    log().end_section(); // end verbose section
  }
  {
    log().subsection().begin<Log::standard>("Standard subsection"); // begin verbose section
    log() << "standard subsection (is printed)" << std::endl;
    log().end_section(); // end verbose section
  }
  log() << "standard section again (is printed)" << std::endl;
  log().end_section(); // end standard section
}

TEST(LogTest, Indentation) {

  // The Log can keep track of indentation levels to create indented sections

  log().set_verbosity(Log::standard);

  log().begin<Log::standard>("Standard section"); // begin standard section
  log() << "standard section (is printed)" << std::endl << std::endl;
  {
    log().increase_indent();
    log().subsection().begin<Log::verbose>("Verbose subsection"); // begin verbose section
    log().indent() << "verbose section (is not printed)" << std::endl;
    log() << std::endl;
    log().end_section(); // end verbose section
    log().decrease_indent();
  }
  {
    log().increase_indent();
    log().subsection().begin<Log::standard>("Standard subsection"); // begin verbose section
    log().indent() << "standard subsection (is printed)" << std::endl;
    log() << std::endl;
    log().end_section(); // end verbose section
    log().decrease_indent();
  }
  log() << "standard section again (is printed)" << std::endl;
  log().end_section(); // end standard section
}

TEST(LogTest, ScopedLogging) {

  // The ScopedLogging class allows swapping what stream CASM::log() and CASM::err_log() refer to
  // for the scope of the ScopedLogging instance

  // write to default log
  log() << "write to default log" << std::endl;

  std::stringstream ss;
  std::stringstream ss_err;
  Log ss_log {ss};
  Log ss_err_log {ss_err};

  {
    ScopedLogging logging {ss_log, ss_err_log};
    log() << "write to stringstream";
    EXPECT_EQ(ss.str(), "write to stringstream");

    err_log() << "write to stringstream";
    EXPECT_EQ(ss_err.str(), "write to stringstream");
  }

  // here, logging has been destructed, so writing to log() no longer appends to ss
  // doesn't get "write to default log again" appended

  log() << "write to default log again" << std::endl;
  EXPECT_EQ(ss.str(), "write to stringstream");

  err_log() << "write to default log again" << std::endl;
  EXPECT_EQ(ss_err.str(), "write to stringstream");

}
