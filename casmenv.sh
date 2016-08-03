# set environment variables recognized by CASM during installation and use

#  Recognized by install scripts and 'casm' CLI executable for CASM install 
#  location (headers and shared libraries)
#  Order of precedence:
#    1) $CASM_PREFIX
#    2) "/usr/local"
#
#export CASM_PREFIX="/usr/local"


#  Recognized by install scripts and 'casm' CLI executable for locating boost 
#  install location
#  Order of precedence:
#    1) $CASM_BOOST_PREFIX
#    2) "" (default uses system defaults)
#
#export CASM_BOOST_PREFIX=""

#  Recognized by install scripts. Use this if linking to boost libraries compiled 
#  without c++11. If defined, (i.e. CASM_BOOST_NO_CXX11_SCOPED_ENUMS=1) will 
#  compile with -DBOOST_NO_CXX11_SCOPED_ENUMS option.
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

#  If "1", sets to compile with debugging symbols. In this case, the optimization 
#  level gets set to -O0, and NDEBUG does not get set. If "0" (default), then 
#  -DNDEBUG is used.
#
#  Order of precedence:
#    1) $CASM_DEBUGSTATE
#    2) "0" (default)
#
#export CASM_DEBUGSTATE="0"

#  Recognized by 'casm' CLI executable for setting runtime compilation
#  order of precedence:
#    1) $CASM_CXXFLAGS
#    2) "-O3 -Wall -fPIC --std=c++11" (default)
#
#export CASM_CXXFLAGS="-O3 -Wall -fPIC --std=c++11"

#  Recognized by 'casm' CLI executable for setting runtime shared library creation options
#  order of precedence:
#    1) $CASM_SOFLAGS
#    2) "-shared -lboost_system" (default)
#
#export CASM_SOFLAGS="-shared -lboost_system"

#  Recognized by 'casm' python module to specify to override 'casm' executable
#  order of precedence:
#    1) $CASM
#    2) 'casm'
#
#export CASM="casm"

#  Recognized by 'casm' python module to specify to 'libcasm' shared library
#  order of precedence:
#    1) $LIBCCASM
#    2) $CASM_PREFIX/lib/libcasm.*'[0]
#    3) '/usr/local/lib/libcasm.*'[0] (default)
#
#export LIBCASM="/usr/local/lib/libcasm.so"

#  Recognized by 'casm' python module to specify to 'libccasm' shared library
#  order of precedence:
#    1) $LIBCCASM
#    2) $CASM_PREFIX/lib/libccasm.*'[0]
#    3) '/usr/local/lib/libccasm.*'[0] (default)
#
#export LIBCCASM="/usr/local/lib/libccasm.so"


#  If CASM_PREFIX is set, make sure our PATH and PYTHONPATH are set
if [ ! -z ${CASM_PREFIX} ]; then
  echo "CASM_PREFIX is set to:"$CASM_PREFIX
  
  # Add $CASM_PREFIX/bin to PATH
  export PATH=$CASM_PREFIX/bin:$PATH
  echo "PATH: "$PATH
  
  # Add $CASM_PREFIX/lib/$PYTHON_VERSION/site-packages to PYTHONPATH
  PYTHON_VERSION=python$(python -c "import sys; print sys.version[:3]")
  export PYTHONPATH=$CASM_PREFIX/lib/$PYTHON_VERSION/site-packages:$PYTHONPATH
  echo "PYTHONPATH: "$PYTHONPATH
  
fi

