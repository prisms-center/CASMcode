#===========================================================#
#To understand why it's structured this way see
#http://karelzak.blogspot.co.uk/2013/02/non-recursive-automake.html

#ACLOCAL_AMFLAGS = -I m4

#===========================================================#
#Start out empty, let the Makemodule.am files fill these in
lib_LTLIBRARIES=
bin_PROGRAMS=
check_PROGRAMS=
man1_MANS=

#===========================================================#
#Define flags for everything ever
include $(srcdir)/common.mk

AM_CXXFLAGS = -DTXT_VERSION='"$(TXT_VERSION)"'\
			  -DEIGEN_DEFAULT_DENSE_INDEX_TYPE=long\
			  -DGZSTREAM_NAMESPACE=gz\
			  -I$(srcdir)/include\
			  -I$(srcdir)/include/casm/external/gzstream\
			  -I$(srcdir)/include/casm/external/qhull/libqhullcpp

AM_CPPFLAGS = -I$(srcdir)/include/casm/external/qhull/libqhull_r/\
			  -I$(srcdir)/include/casm/external/gzstream

#AM_CXXFLAGS+=$(BOOST_CXXFLAGS)? AM_LDFLAGS+=$(BOOST_LDFLAGS)?

#===========================================================#
#Create shared casm libraries
include $(srcdir)/src/casm/Makemodule.am
include $(srcdir)/src/ccasm/Makemodule.am

#===========================================================#
#Take care of all the include headers
include $(srcdir)/include/casm/Makemodule.am
include $(srcdir)/include/ccasm/Makemodule.am

#===========================================================#
#Take care of applications
include $(srcdir)/apps/casm/Makemodule.am

#===========================================================#
#Take care of tests
include $(srcdir)/tests/unit/Makemodule.am

#===========================================================#
#Misc targets
.FORCE:

#Files that arent sources for executalbes or headers that get installed, but are needed in the distribution.
EXTRA_DIST=./build-aux/casm_version.txt\
		   ./man/casm.1\
		   ./python/Makefile.in\
		   ./python/Makefile.am\
		   ./tests/unit/Common.hh

