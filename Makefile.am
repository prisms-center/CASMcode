#===========================================================#
#To understand why it's structured this way see
#http://karelzak.blogspot.co.uk/2013/02/non-recursive-automake.html

ACLOCAL_AMFLAGS = -I build-aux/m4

#===========================================================#
#Start out empty, let the Makemodule.am files fill these in
lib_LTLIBRARIES=
bin_PROGRAMS=
check_PROGRAMS=
man1_MANS=

#===========================================================#
#Define flags for everything ever
TXT_VERSION:=$(shell cat $(top_srcdir)/build-aux/casm_version.txt)

AM_CXXFLAGS = -DTXT_VERSION='"$(TXT_VERSION)"'\
			  -DEIGEN_DEFAULT_DENSE_INDEX_TYPE=long\
			  -DGZSTREAM_NAMESPACE=gz\
			  -I$(srcdir)/include\
			  -I$(srcdir)/include/casm/external/gzstream\
			  -I$(srcdir)/include/casm/external/qhull/libqhullcpp

AM_CPPFLAGS = -I$(srcdir)/include/casm/external/qhull/libqhull_r/\
			  -I$(srcdir)/include/casm/external/gzstream\
			  $(BOOST_CPPFLAGS)

AM_LDFLAGS = $(BOOST_LDFLAGS)

#AM_CXXFLAGS+=$(BOOST_CXXFLAGS)? AM_LDFLAGS+=$(BOOST_LDFLAGS)?

#============================================================#
#Always have included modules use paths relative to this
#parent Makefile.am. Avoid variable expansion in *_SOURCES
#============================================================#
#Create shared casm libraries
include $(srcdir)/src/casm/Makemodule.am
include $(srcdir)/src/ccasm/Makemodule.am

#============================================================#
#Take care of all the include headers
include $(srcdir)/include/casm/Makemodule.am
include $(srcdir)/include/ccasm/Makemodule.am

#============================================================#
#Take care of applications
include $(srcdir)/apps/casm/Makemodule.am

#============================================================#
#Take care of tests
include $(srcdir)/tests/unit/Makemodule.am

#============================================================#
#Misc targets
.FORCE:

test:
	echo $(srcdir)
	echo $(top_srcdir)
	echo $(AM_LDFLAGS)
	echo $(AM_CXXFLAGS)
	echo $(AM_CPPFLAGS)
	echo $(casm_LDFLAGS)
	echo $(BOOST_LDFLAGS)
	echo $(BOOST_CXXFLAGS)
	echo $(BOOST_CPPFLAGS)
	echo $(TXT_VERSION)

#============================================================#
#Recurse into these places (must be direct children, don't be fancy)
#SUBDIRS= python

#============================================================#
#Files that arent sources for executalbes or headers that get installed, but are needed in the distribution.
#Even though the include/casm/external is distributed through headers, the licensing and other files should
#also be included.
EXTRA_DIST=./build-aux/casm_version.txt\
		   ./man\
		   ./python\
		   ./tests\
		   ./include/casm/external

#Extra files that should be cleaned up upon `distclean`
DISTCLEANFILES=./python/Makefile

