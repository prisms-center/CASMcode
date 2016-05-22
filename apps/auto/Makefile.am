AM_CXXFLAGS = $(BOOST_CXXFLAGS) -Wno-deprecated-register -Wno-deprecated-declarations -DEIGEN_DEFAULT_DENSE_INDEX_TYPE=long -Wno-unused-parameter -DGZSTREAM_NAMESPACE=gz -I../../include
AM_LDFLAGS = $(BOOST_LDFLAGS)
bin_PROGRAMS = autotest
autotest_LDADD=../../src/casm/libcasm.la $(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_REGEX_LIB)

autotest_SOURCES = auto.cpp
