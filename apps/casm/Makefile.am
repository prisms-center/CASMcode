AM_CXXFLAGS = $(BOOST_CXXFLAGS) -DEIGEN_DEFAULT_DENSE_INDEX_TYPE=long -DGZSTREAM_NAMESPACE=gz -I$(top_srcdir)/include
AM_LDFLAGS = $(BOOST_LDFLAGS)

bin_PROGRAMS = casm
casm_SOURCES = casm.cpp bset.cc composition.cc enum.cc format.cc import.cc init.cc monte.cc perturb.cc ref.cc run.cc select.cc settings.cc status.cc super.cc sym.cc update.cc view.cc 
casm_LDADD=$(top_srcdir)/src/casm/libcasm.la $(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_REGEX_LIB) $(BOOST_CHRONO_LIB)

man1_MANS = $(top_srcdir)/man/casm.1
