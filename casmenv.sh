### Set environment variables recognized by CASM during installation and use

#  Recognized by scons install scripts and 'casm' CLI executable for CASM install
#  location assuming $CASM_PREFIX/include for headers and $CASM_PREFIX/lib for libcasm
#  Order of precedence:
#    1) $CASM_INCLUDEDIR, $CASM_LIBDIR
#    2) $CASM_PREFIX/include, $CASM_PREFIX/lib
#    3) "/usr/local"
#
#export CASM_PREFIX="/usr/local"
#export CASM_INCLUDEDIR="/usr/local/include"
#export CASM_LIBDIR="/usr/local/lib"


#  Recognized by install scripts and 'casm' CLI executable for locating boost 
#  install location
#  Order of precedence:
#    1) $CASM_BOOST_INCLUDEDIR, $CASM_BOOST_LIBDIR
#    2) $CASM_BOOST_PREFIX/include $CASM_BOOST_PREFIX/lib
#    3) "" (default uses system defaults)
#
#export CASM_BOOST_PREFIX=""
#export CASM_BOOST_INCLUDEDIR=""
#export CASM_BOOST_LIBDIR=""
 

#  Recognized by install scripts. Use this if linking to boost libraries compiled without c++11. If defined, (i.e. CASM_BOOST_NO_CXX11_SCOPED_ENUMS=1) will compile with -DBOOST_NO_CXX11_SCOPED_ENUMS option.
#  Order of precedence:
#    1) if $CASM_BOOST_NO_CXX11_SCOPED_ENUMS defined
#    2) not defined (default)
#
#export CASM_BOOST_NO_CXX11_SCOPED_ENUMS="1"

#  Recognized by install scripts and 'casm' CLI executable for setting compiler
#  order of precedence:
#    1) $CASM_CXX
#    2) $CXX
#    3) "g++"
#
#export CASM_CXX="g++"

#  Sets the -O optimization compiler option. Default is '3', which indicates -O3.
#  Order of precedence:
#    1) $CASM_OPTIMIZATIONLEVEL
#    2) "3" (default)
#
#export CASM_OPTIMIZATIONLEVEL="3"

#   If "1", sets to compile with debugging symbols. In this case, the optimization level gets set to -O0, and NDEBUG does not get set. If "0" (default), then -DNDEBUG is used.
#
#  Order of precedence:
#    1) $CASM_DEBUGSTATE
#    2) "0" (default)
#
#export CASM_DEBUGSTATE="0"

#  recognized by 'casm' CLI executable for setting runtime compilation
#  order of precedence:
#    1) $CASM_CXXFLAGS
#    2) "-O3 -Wall -fPIC --std=c++11" (default)
#
#export CASM_CXXFLAGS="-O3 -Wall -fPIC --std=c++11"

#  recognized by 'casm' CLI executable for setting runtime shared library creation options
#  order of precedence:
#    1) $CASM_SOFLAGS
#    2) "-shared -lboost_system" (default)
#
#export CASM_SOFLAGS="-shared -lboost_system"

#  recognized by 'casm' python module to specify to override 'casm' executable
#  order of precedence:
#    1) $CASM
#    2) 'casm'
#
#export CASM="casm"

#  recognized by 'casm' python module to specify to 'libcasm' shared library
#  order of precedence:
#    1) $LIBCCASM
#    2) $CASM_PREFIX/lib/libcasm.*'[0]
#    3) '/usr/local/lib/libcasm.*'[0] (default)
#
#export LIBCASM="/usr/local/lib/libcasm.so"

#  recognized by 'casm' python module to specify to 'libccasm' shared library
#  order of precedence:
#    1) $LIBCCASM
#    2) $CASM_PREFIX/lib/libccasm.*'[0]
#    3) '/usr/local/lib/libccasm.*'[0] (default)
#
#export LIBCCASM="/usr/local/lib/libccasm.so"


### Testing environment

#  To run the non-installed casm for testing purpose, set this to the location
#  of the CASM git repository 
#
#CASM_REPO="/path/to/CASMcode"



#### You probably don't need to change these ####

#  If CASM_PREFIX is set, update paths
if [ ! -z ${CASM_PREFIX} ]; then
  # Add $CASM_PREFIX/bin to PATH
  export PATH=$CASM_PREFIX/bin:$PATH
  
  # Add $CASM_PREFIX/lib/$PYTHON_VERSION/site-packages to PYTHONPATH
  PYTHON_VERSION=python$(python -c "import sys; print sys.version[:3]")
  export PYTHONPATH=$CASM_PREFIX/lib/$PYTHON_VERSION/site-packages:$PYTHONPATH

fi

#  If CASM_BOOST_PREFIX is set, update library search path
if [ ! -z ${CASM_BOOST_PREFIX} ]; then
  
  # For Linux, set LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=$CASM_BOOST_PREFIX/lib:$LD_LIBRARY_PATH
  
  # For Mac, set DYLD_LIBRARY_FALLBACK_PATH
  export DYLD_FALLBACK_LIBRARY_PATH=$CASM_BOOST_PREFIX/lib:$DYLD_FALLBACK_LIBRARY_PATH
  
fi

# If testing:
if [ ! -z ${CASM_REPO} ]; then
  
  export CASM=$CASM_REPO/bin/casm
  export LIBCASM=$CASM_REPO/lib/libcasm.dylib
  export LIBCCASM=$CASM_REPO/lib/libccasm.dylib
  export PATH=$CASM_REPO/bin:$CASM_REPO/python/casm/scripts:$PATH
  export PYTHONPATH=$CASM_REPO/python/casm:$PYTHONPATH
  
  # For testing on Linux, use LD_LIBRARY_PATH:
  export LD_LIBRARY_PATH=$CASM_REPO/lib:$LD_LIBRARY_PATH
  
  # For testing on Mac, use DYLD_FALLBACK_LIBRARY_PATH:
  export DYLD_FALLBACK_LIBRARY_PATH=$CASM_REPO/lib:$DYLD_FALLBACK_LIBRARY_PATH
  
fi


