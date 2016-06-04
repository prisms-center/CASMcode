AM_CXXFLAGS = $(BOOST_CXXFLAGS) -I$(top_srcdir)/include
AM_LDFLAGS = $(top_srcdir)/src/casm/libcasm.la $(BOOST_LDFLAGS) 

check_PROGRAMS = casm_units
casm_units_SOURCES = ./clex/test_Clexulator.cc ./Common.cpp ./unit_test.cpp ./App/App_test.cpp ./system/RuntimeLibrary_test.cpp ./casm_io/jsonParser_test.cpp ./container/Counter_test.cpp ./clex/ConfigIO_test.cpp ./clex/ChemicalReference_test.cpp ./clex/NeighborList_test.cpp ./clex/PrimClex_test.cpp ./clex/Clexulator_test.cpp ./crystallography/BasicStructure_test.cpp ./crystallography/Lattice_test.cpp ./crystallography/Coordinate_test.cpp ./crystallography/Structure_test.cpp 
casm_units_LDADD=$(top_srcdir)/src/casm/libcasm.la $(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_REGEX_LIB) $(BOOST_CHRONO_LIB) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) $(DL_LIB)
