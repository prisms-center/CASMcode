# http://www.scons.org/doc/production/HTML/scons-user.html
# This is: Sconstruct

import os, glob, copy, shutil, subprocess

from os.path import join

Help("""
      Type: 'scons' to build all binaries,
            'scons install' to install all libraries, binaries, scripts and python packages,
            'scons test' to run all tests,
            'scons unit' to run all unit tests,
            'scons A_UNIT_TEST' to run a particular unit test (where A_UNIT_TEST 
                                is replaced with the name of the particular unit test, 
                                typically a class name),
            'scons casm_test' to run tests/casm.
            
      In all cases, add '-c' to perform a clean up or uninstall.
      
      Default compile options are: '-O3 -DNDEBUG -Wno-unused-parameter'
      
      
      Recognized environment variables:
      
      $CXX:
        Explicitly set the C++ compiler. If not set, scons chooses a default compiler.
      
      $CASMPREFIX:
        Where to install CASM. By default, this uses '/usr/local'. Then header files are
        installed in '$CASMPREFIX/include', shared libraries in '$CASMPREFIX/lib', executables
        in '$CASMPREFIX/bin', and the path is used for the setup.py --prefix option for 
        installing python packages.
      
      $CASMBOOST_PATH:
        Search path for Boost. '$CASMBOOST_PATH/include' is searched for header files, and
        '$CASMBOOST_PATH/lib' for libraries. Boost and CASM should be compiled with the 
        same compiler.

      $OPTIMIZATIONLEVEL:
        Sets the -O optimization compiler option. If not set, uses -O3.

      $DEBUGSTATE:
        Sets to compile with debugging symbols. In this case, the optimization level gets 
        set to -O0, and NDEBUG does not get set.

      $LD_LIBRARY_PATH:
        Search path for dynamic libraries, may need $CASMBOOST_LIB 
        and $CASMPREFIX/lib added to it.
        On Mac OS X, this variable is $DYLD_FALLBACK_LIBRARY_PATH.
        This should be added to your ~/.bash_profile (Linux) or ~/.profile (Mac).
      
      $CASMBOOST_NO_CXX11_SCOPED_ENUMS:
        If defined, will compile with -DCASMBOOST_NO_CXX11_SCOPED_ENUMS. Use this
        if linking to boost libraries compiled without c++11.
      
      
      Additional options that override environment variables:
      
      Use 'cxx=X' to set the C++ compiler. Default is chosen by scons.
          'opt=X' to set optimization level, '-OX'. Default is 3.
          'debug=X' with X=0 to use '-DNDEBUG', 
             or with X=1 to set debug mode compiler options '-O0 -g -save-temps'.
             Overrides $DEBUGSTATE.
          'prefix=X' to set installation directory. Default is '/usr/local'. Overrides $CASMPREFIX.
          'boost_path=X' set boost search path. Overrides $CASMBOOST_PATH.
          'boost_no_cxx11_scoped_enums=1' to use '-DBOOST_NO_CXX11_SCOPED_ENUMS'.
             Overrides $CASMBOOST_NO_CXX11_SCOPED_ENUMS.
     """)

def version(version_number):
  
  # check if git installed
  try:
    # pipe output to /dev/null for silence
    null = open("/dev/null", "w")
    subprocess.Popen("git", stdout=null, stderr=null)
    null.close()

  except OSError:
    return version_number

  # want to get the current git branch name, if in a git repository, else ''
  process = subprocess.Popen(['git', 'rev-parse', '--abbrev-ref', 'HEAD'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  
  # get the stdout, split if any '/', and use the last bit (the branch name), strip any extra space
  branch = process.communicate()[0].split('/')[-1].strip()

  if branch == '':
    return version_number

  # when compiling from a git repo use a developement version number
  # which contains the branch name, short hash, and tag (if tagged)
    
  # get the short hash
  process = subprocess.Popen('git rev-parse --short HEAD'.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  commit = process.communicate()[0].strip()
  version_number += '+' + commit
  
  # check if tracked files have changes, if so, add '+changes'
  process = subprocess.Popen('git status --porcelain --untracked-files=no'.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  changes = process.communicate()[0].strip()
  if changes != '':
    version_number += ".changes"
  
  return version_number
  


##### Set version_number

version_number = version('0.2a0')
url = 'https://github.com/prisms-center/CASMcode'
Export('version_number', 'url')


##### Environment setup

# include paths
include_paths = [join(os.getcwd(), i) for i in ['include', 'include/casm/app'] ]

# lib paths
lib_paths = []

# command-line variables (C and C++)
ccflags = []
cxxflags = []

#ccflags.append('-Wall')

ccflags.append('-Wno-unused-parameter')

debug_level = '0'
if 'debug' in ARGUMENTS:
  debug_level = ARGUMENTS.get('debug')
elif 'DEBUGSTATE' in os.environ:
  debug_level = os.environ['DEBUGSTATE']

if debug_level == '0':
  ccflags = ccflags + ['-DNDEBUG']
  cxxflags = cxxflags + ['-DNDEBUG']
elif debug_level == '1':
  ccflags = ccflags + ['-g', '-save-temps']
  cxxflags = cxxflags + ['-g', '-save-temps']  

if 'OPTIMIZATIONLEVEL' in os.environ:
  opt_level = os.environ['OPTIMIZATIONLEVEL']
else:
  if debug_level == '0':
    opt_level = '3'
  else:
    opt_level = '0'

ccflags.append("-O" + ARGUMENTS.get('opt',opt_level))
cxxflags.append("-O" + ARGUMENTS.get('opt',opt_level))

# C++ only
#cxxflags = []
cxxflags.append('--std=c++11')
cxxflags.append('-Wno-deprecated-register')
cxxflags.append('-Wno-deprecated-declarations')
cxxflags.append('-DEIGEN_DEFAULT_DENSE_INDEX_TYPE=long')

boost_no_cxx11_scoped_enums = ARGUMENTS.get('boost_no_cxx11_scoped_enums', '0')
if 'CASMBOOST_NO_CXX11_SCOPED_ENUMS' in os.environ:
  boost_no_cxx11_scoped_enums = '1'
if boost_no_cxx11_scoped_enums == '1':
  cxxflags.append('-DBOOST_NO_CXX11_SCOPED_ENUMS')

# set gzstream namespace to 'gz'
ccflags.append('-DGZSTREAM_NAMESPACE=gz')


boost_path = None
if 'boost_path' in ARGUMENTS:
  boost_path = ARGUMENTS.get('boost_path')
elif 'CASMBOOST_PATH' in os.environ:
  boost_path = os.environ['CASMBOOST_PATH']
if(boost_path != None):
  include_paths.append(os.path.join(boost_path, 'include'))
  lib_paths.append(os.path.join(boost_path, 'lib'))

# where everything is built
build_lib_paths = copy.deepcopy(lib_paths)
build_lib_paths.append(os.path.join(os.getcwd(), 'lib'))
Export('build_lib_paths')


# where everything should be installed
install_lib_paths = copy.deepcopy(lib_paths)
prefix = '/usr/local'
if 'prefix' in ARGUMENTS:
  prefix = ARGUMENTS.get('prefix')
elif 'CASMPREFIX' in os.environ:
  prefix = os.environ['CASMPREFIX']
install_lib_paths.append(os.path.join(prefix, 'lib'))
Export('install_lib_paths')

env = Environment(ENV = os.environ,
                  CCFLAGS = ccflags,
                  CXXFLAGS = cxxflags,
                  CPPPATH = include_paths,
                  LIBPATH = lib_paths,
                  PREFIX = prefix)

# set a non-default c++ compiler
if 'CXX' in os.environ:
  env.Replace(CXX = env['CXX'])
elif 'cxx' in ARGUMENTS:
  env.Replace(CXX = ARGUMENTS.get('cxx'))

# where the shared libraries should go
env.Append(CASM_LIB = os.path.join(os.getcwd(), 'lib'))

# where the compiled binary should go
unit_test_bin = os.path.join(os.getcwd(), 'tests', 'unit', 'bin')
env.Append(UNIT_TEST_BIN = unit_test_bin)
casm_bin = os.path.join(os.getcwd(), 'bin')
env.Append(CASM_BIN = casm_bin)

# collect header files
env.Append(CASM_SOBJ = [])

# collect everything that will go into the casm library
env.Append(CASM_SOBJ = [])

# collect everything that will go into the casm c library
env.Append(CCASM_SOBJ = [])

# Whenever a new Alias is declared, provide a check for if a 'test' or 'installation' is being done
# and store the result in these environment variables, then use this to prevent undesired clean up
env.Append(COMPILE_TARGETS = [])
env.Append(INSTALL_TARGETS = [])
env.Append(IS_TEST = 0)
env.Append(IS_INSTALL = 0)

# make compiler errors and warnings in color
env['ENV']['TERM'] = os.environ['TERM']

# set testing environment
env['ENV']['PATH'] = env['CASM_BIN'] + ":" + env['ENV']['PATH']


##### Call all SConscript files for shared objects

# build src/casm/external/gzstream
SConscript(['src/casm/external/gzstream/SConscript'], {'env':env})

# build src/casm
SConscript(['src/casm/SConscript'], {'env':env})

# build src/ccasm
SConscript(['src/ccasm/SConscript'], {'env':env})


##### Make single dynamic library 

linkflags = ""
if env['PLATFORM'] == 'darwin':
  linkflags = ['-install_name', '@rpath/libcasm.dylib']

# use boost libraries
boost_libs = ['boost_system', 'boost_filesystem', 'boost_program_options', 'boost_regex', 'boost_chrono']

# build casm shared library from all shared objects
casm_lib = env.SharedLibrary(os.path.join(env['CASM_LIB'], 'casm'), 
                             env['CASM_SOBJ'], 
                             LIBPATH=build_lib_paths,
                             LINKFLAGS=linkflags,
                             LIBS=boost_libs + ['z'])
                             
env['COMPILE_TARGETS'] = env['COMPILE_TARGETS'] + casm_lib
Export('casm_lib')
env.Alias('libcasm', casm_lib)

# Library Install instructions
casm_lib_install = env.SharedLibrary(os.path.join(env['PREFIX'], 'lib', 'casm'), 
                                     env['CASM_SOBJ'], 
                                     LIBPATH=install_lib_paths, 
                                     LINKFLAGS=linkflags,
                                     LIBS=boost_libs + ['z'])
Export('casm_lib_install')
env.Alias('casm_lib_install', casm_lib_install)
env['INSTALL_TARGETS'] = env['INSTALL_TARGETS'] + [casm_lib_install]

if 'casm_lib_install' in COMMAND_LINE_TARGETS:
    env['IS_INSTALL'] = 1


##### extern C dynamic library 

linkflags = ""
if env['PLATFORM'] == 'darwin':
  linkflags = ['-install_name', '@rpath/libccasm.dylib']

# build casm shared library from all shared objects
ccasm_lib = env.SharedLibrary(os.path.join(env['CASM_LIB'], 'ccasm'), 
                             env['CCASM_SOBJ'], 
                             LIBPATH=build_lib_paths,
                             LINKFLAGS=linkflags,
                             LIBS=boost_libs + ['casm'])
                             
env['COMPILE_TARGETS'] = env['COMPILE_TARGETS'] + ccasm_lib
Export('ccasm_lib')
env.Alias('libccasm', ccasm_lib)

# Library Install instructions
ccasm_lib_install = env.SharedLibrary(os.path.join(env['PREFIX'], 'lib', 'ccasm'), 
                                     env['CCASM_SOBJ'], 
                                     LIBPATH=install_lib_paths, 
                                     LINKFLAGS=linkflags,
                                     LIBS=boost_libs + ['casm'])
Export('ccasm_lib_install')
env.Alias('ccasm_lib_install', ccasm_lib_install)
env['INSTALL_TARGETS'] = env['INSTALL_TARGETS'] + [ccasm_lib_install]

if 'ccasm_lib_install' in COMMAND_LINE_TARGETS:
    env['IS_INSTALL'] = 1


# Include Install instructions

# for all header files: install h/path/to/X.hh as prefix/casm/path/to/X.hh
# we'll create a dict of file -> install dir

include_dir = os.path.join(env['PREFIX'],'include')
casm_include = os.path.join('include', 'casm')
Export('casm_include')

casm_include_install = env.Install(include_dir, casm_include)

Export('casm_include_install')
env.Alias('casm_include_install', casm_include_install)
env.Clean('casm_include_install', os.path.join(include_dir,'casm'))
env['INSTALL_TARGETS'] = env['INSTALL_TARGETS'] + casm_include_install

if 'casm_include_install' in COMMAND_LINE_TARGETS:
  env['IS_INSTALL'] = 1

##### Call all SConscript files for executables

# build apps/casm
SConscript(['apps/casm/SConscript'], {'env':env})

# tests/unit
SConscript(['tests/unit/SConscript'], {'env': env})

# tests/casm
SConscript(['tests/casm/SConscript'], {'env': env})


##### Python packages

# install python packages and scripts
SConscript(['python/casm/SConscript'], {'env':env})
SConscript(['python/vasp/SConscript'], {'env':env})



##### Make combined alias 'test'

# Execute 'scons test' to compile & run integration and unit tests
env.Alias('test', ['unit', 'casm_test'])

if 'test' in COMMAND_LINE_TARGETS:
    env['IS_TEST'] = 1


##### Make combined alias 'install'

# Execute 'scons install' to install all binaries, scripts and python modules
installable = ['casm_include_install', 'casm_lib_install', 'ccasm_lib_install', 'casm_install', 'pycasm_install', 'pyvasp_install']
env.Alias('install', installable)

if 'install' in COMMAND_LINE_TARGETS:
    env['IS_INSTALL'] = 1

# scons is not checking if any header files changed
# if we're supposed to install them, rm -r include_dir/casm
would_install_include = ['casm_include_install', 'install']
if [i for i in would_install_include if i in COMMAND_LINE_TARGETS]:
  path = os.path.join(include_dir, 'casm')
  if os.path.exists(path):
    print "rm", path
    shutil.rmtree(path)

##### Clean up instructions

if env['IS_INSTALL'] or env['IS_TEST']:
  env.NoClean(env['COMPILE_TARGETS'])
if not env['IS_INSTALL']:
  env.NoClean(env['INSTALL_TARGETS'])
  if debug_level == '1':
    env.Clean(casm_lib, Glob('*.s') + Glob('*.ii'))
