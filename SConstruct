# http://www.scons.org/doc/production/HTML/scons-user.html
# This is: Sconstruct

import sys, os, glob, copy, shutil, subprocess, imp, re

from os.path import join

Help("""
      Type: 'scons configure' to run configuration checks,
            'scons' to build all binaries,
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
      
      $CASM_CXX, $CXX:
        Explicitly set the C++ compiler. If not set, scons chooses a default compiler.
      
      $CASM_PREFIX:
        Where to install CASM. By default, this uses '/usr/local'. Then header files are
        installed in '$CASM_PREFIX/include', shared libraries in '$CASM_PREFIX/lib', executables
        in '$CASM_PREFIX/bin', and the path is used for the setup.py --prefix option for 
        installing python packages.
      
      $CASM_BOOST_PREFIX:
        Search path for Boost. '$CASM_BOOST_PREFIX/include' is searched for header files, and
        '$CASM_BOOST_PREFIX/lib' for libraries. Boost and CASM should be compiled with the 
        same compiler.

      $CASM_OPTIMIZATIONLEVEL:
        Sets the -O optimization compiler option. If not set, uses -O3.

      $CASM_DEBUGSTATE:
        Sets to compile with debugging symbols. In this case, the optimization level gets 
        set to -O0, and NDEBUG does not get set.

      $LD_LIBRARY_PATH (Linux) or $DYLD_FALLBACK_LIBRARY_PATH (Mac):
        Search path for dynamic libraries, may need $CASM_BOOST_PREFIX/lib 
        and $CASM_PREFIX/lib added to it.
        This should be added to your ~/.bash_profile (Linux) or ~/.profile (Mac).
      
      $CASM_BOOST_NO_CXX11_SCOPED_ENUMS:
        If defined, will compile with -DCASM_BOOST_NO_CXX11_SCOPED_ENUMS. Use this
        if linking to boost libraries compiled without c++11.
      
      
      Additional options that override environment variables:
      
      Use 'cxx=X' to set the C++ compiler. Default is chosen by scons.
          'opt=X' to set optimization level, '-OX'. Default is 3.
          'debug=X' with X=0 to use '-DNDEBUG', 
             or with X=1 to set debug mode compiler options '-O0 -g -save-temps'.
             Overrides $CASM_DEBUGSTATE.
          'prefix=X' to set installation directory. Default is '/usr/local'. Overrides $CASM_PREFIX.
          'boost_prefix=X' set boost search path. Overrides $CASM_BOOST_PPREFIX.
          'boost_no_cxx11_scoped_enums=1' to use '-DBOOST_NO_CXX11_SCOPED_ENUMS'.
             Overrides $CASM_BOOST_NO_CXX11_SCOPED_ENUMS.
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

def cxx():
  """Get c++ compiler"""
  # set a non-default c++ compiler
  if 'cxx' in ARGUMENTS:
    return ARGUMENTS.get('cxx')
  elif 'CASM_CXX' in os.environ:
    return os.environ['CASM_CXX']
  elif 'CXX' in os.environ:
    return os.environ['CXX']
  return "g++"

def prefix():
  """Install location"""
  if 'prefix' in ARGUMENTS:
    return ARGUMENTS.get('prefix')
  elif 'CASM_PREFIX' in os.environ:
    return os.environ['CASM_PREFIX']
  return  '/usr/local'

def boost_prefix():
  """Boost location"""
  if 'boost_prefix' in ARGUMENTS:
    return ARGUMENTS.get('boost_prefix')
  elif 'CASM_BOOST_PREFIX' in os.environ:
    return os.environ['CASM_BOOST_PREFIX']
  return None

def include_path(self, prefix, name):
  """
  Return prefix/include where prefix/include/name exists, or None
  
  If prefix is None, check /usr and /usr/local, else just check prefix
  
  """
  if prefix is None:
    for path in ['/usr', '/usr/local']:
      dirpath = join(path, 'include', name)
      if os.path.isdir(dirpath) and os.access(dirpath, os.R_OK):
        return join(path, 'include')
  else:
    dirpath = join(prefix, 'include', name)
    if os.path.isdir(dirpath) and os.access(dirpath, os.R_OK):
      return join(prefix, 'include')
  return None

def lib_path(prefix, name):
  """
  Return prefix/lib where prefix/include/name exists, or None
  
  If prefix is None, check /usr and /usr/local, else just check prefix
  
  """
  if prefix is None:
    for path in ['/usr', '/usr/local']:
      dirpath = join(path, 'include', name)
      if os.path.isdir(dirpath) and os.access(dirpath, os.R_OK):
        return join(path, 'lib')
  else:
    dirpath = join(prefix, 'include', name)
    if os.path.isdir(dirpath) and os.access(dirpath, os.R_OK):
      return join(prefix, 'lib')
  return None

def lib_name(prefix, includename, libname):
  """
  Uses re.match('lib(' + libname + '.*)\.(dylib|a|so).*',string) on all files
  in the prefix/lib directory to get the libname to use. If none found, return None.
  """
  for p in os.listdir(lib_path(prefix, includename)):
    m = re.match('lib(' + libname + '.*)\.(dylib|a|so).*',p)
    if m:
      return m.group(1)
  return None

def debug_level():
  if 'debug' in ARGUMENTS:
    return ARGUMENTS.get('debug')
  elif 'CASM_DEBUGSTATE' in os.environ:
    return os.environ['CASM_DEBUGSTATE']
  return '0'

def opt_level():
  if 'opt' in ARGUMENTS:
    return ARGUMENTS.get('opt')
  elif 'CASM_OPTIMIZATIONLEVEL' in os.environ:
    return os.environ['CASM_OPTIMIZATIONLEVEL']
  else:
    _debug_level = debug_level()
    if _debug_level == '0':
      return '3'
    else:
      return '0'

def boost_no_cxx11_scoped_enums():
  if 'boost_no_cxx11_scoped_enums' in ARGUMENTS:
    return True
  elif 'CASM_BOOST_NO_CXX11_SCOPED_ENUMS' in os.environ:
    return True
  return False

def compile_flags():
  # command-line variables (C and C++)
  ccflags = []
  cxxflags = []

  #ccflags.append('-Wall')
  ccflags.append('-Wno-unused-parameter')

  _debug_level = debug_level()
  
  if _debug_level == '0':
    ccflags = ccflags + ['-DNDEBUG']
    cxxflags = cxxflags + ['-DNDEBUG']
  elif _debug_level == '1':
    ccflags = ccflags + ['-g', '-save-temps']
    cxxflags = cxxflags + ['-g', '-save-temps']  

  _opt = "-O" + opt_level()

  ccflags.append(_opt)
  cxxflags.append(_opt)

  # C++ only
  #cxxflags = []
  cxxflags.append('--std=c++11')
  cxxflags.append('-Wno-deprecated-register')
  cxxflags.append('-Wno-deprecated-declarations')
  cxxflags.append('-DEIGEN_DEFAULT_DENSE_INDEX_TYPE=long')

  if boost_no_cxx11_scoped_enums():
    cxxflags.append('-DBOOST_NO_CXX11_SCOPED_ENUMS')
  
  # set gzstream namespace to 'gz'
  ccflags.append('-DGZSTREAM_NAMESPACE=gz')
  
  return (ccflags, cxxflags)

##### Set version_number

version_number = version('0.3a0')
url = 'https://github.com/prisms-center/CASMcode'
Export('version_number', 'url')


##### Environment setup

env = Environment()

# C++ compiler
env.Replace(CXX=cxx())

# C and C++ flags
ccflags, cxxflags = compile_flags()
env.Replace(CCFLAGS=ccflags, CXXFLAGS=cxxflags)

# where the shared libraries should go
env.Append(LIBDIR = join(os.getcwd(), 'lib'))

# where the compiled binary should go
env.Append(BINDIR = join(os.getcwd(), 'bin'))

# where the headers are
env.Append(INCDIR = join(os.getcwd(), 'include'))

# where test binaries should go
env.Append(UNIT_TEST_BINDIR = join(os.getcwd(), 'tests', 'unit', 'bin'))

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

# collect casm prefix
env.Replace(PREFIX=prefix())

# collect boost prefix
env.Replace(BOOST_PREFIX=boost_prefix())

# collect library names
env.Replace(z='z') 
env.Replace(dl='dl') 

boost_libs = [
  'boost_system', 
  'boost_filesystem', 
  'boost_program_options',
  'boost_regex',
  'boost_chrono',
  'boost_timer',
  'boost_unit_test_framework'
]
for x in boost_libs:
  env[x] = lib_name(env['BOOST_PREFIX'], 'boost', x)

# make compiler errors and warnings in color
env['ENV']['TERM'] = os.environ['TERM']

# set testing environment (for running tests)
env['ENV']['PATH'] = env['BINDIR'] + ":" + env['ENV']['PATH']

# set LD_LIBRARY_PATH or DYLD_FALLBACK_LIBRARY_PATH (for running configuration tests)
for x in ['LD_LIBRARY_PATH', 'DYLD_LIBRARY_FALLBACK_PATH']:
  if x in os.environ:
    env['ENV'][x] = os.environ[x]

# add methods to use elsewhre
env.AddMethod(include_path)

##### Call all SConscript files for shared objects

# build src/casm/external/gzstream
SConscript(['src/casm/external/gzstream/SConscript'], {'env':env})

# build src/casm
SConscript(['src/casm/SConscript'], {'env':env})

# build src/ccasm
SConscript(['src/ccasm/SConscript'], {'env':env})


##### Make single dynamic library 

# include paths
candidates = [
  env['INCDIR'], 
  env.include_path(env['BOOST_PREFIX'], 'boost')
]
cpppath = [x for x in candidates if x is not None]

# link paths
build_lib_paths = [env['LIBDIR']]
if 'BOOST_PREFIX' in env and env['BOOST_PREFIX'] is not None:
  build_lib_paths.append(join(env['BOOST_PREFIX'], 'lib'))
Export('build_lib_paths')

# link flags
linkflags = ""
if env['PLATFORM'] == 'darwin':
  linkflags = ['-install_name', '@rpath/libcasm.dylib']

# use boost libraries
libs = [
  env['boost_system'], 
  env['boost_filesystem'], 
  env['boost_program_options'], 
  env['boost_regex'], 
  env['boost_chrono'],
  env['z']]

# build casm shared library from all shared objects
casm_lib = env.SharedLibrary(join(env['LIBDIR'], 'casm'), 
                             env['CASM_SOBJ'], 
                             CPPPATH=cpppath,
                             LIBPATH=build_lib_paths,
                             LINKFLAGS=linkflags,
                             LIBS=libs)
                             
env['COMPILE_TARGETS'] = env['COMPILE_TARGETS'] + casm_lib
Export('casm_lib')
env.Alias('libcasm', casm_lib)

# Library Install instructions
casm_lib_install = env.SharedLibrary(join(env['PREFIX'], 'lib', 'casm'), 
                                     env['CASM_SOBJ'], 
                                     CPPPATH=cpppath,
                                     LIBPATH=build_lib_paths, 
                                     LINKFLAGS=linkflags,
                                     LIBS=libs)
Export('casm_lib_install')
env.Alias('casm_lib_install', casm_lib_install)
env['INSTALL_TARGETS'] = env['INSTALL_TARGETS'] + [casm_lib_install]

if 'casm_lib_install' in COMMAND_LINE_TARGETS:
    env['IS_INSTALL'] = 1


##### extern C dynamic library 

linkflags = ""
if env['PLATFORM'] == 'darwin':
  linkflags = ['-install_name', '@rpath/libccasm.dylib']

install_lib_paths = build_lib_paths + [join(env['PREFIX'], 'lib')]

# build casm shared library from all shared objects
ccasm_lib = env.SharedLibrary(join(env['LIBDIR'], 'ccasm'), 
                             env['CCASM_SOBJ'], 
                             CPPPATH=cpppath,
                             LIBPATH=install_lib_paths,
                             LINKFLAGS=linkflags,
                             LIBS=libs + ['casm'])
                             
env['COMPILE_TARGETS'] = env['COMPILE_TARGETS'] + ccasm_lib
Export('ccasm_lib')
env.Alias('libccasm', ccasm_lib)

# Library Install instructions
ccasm_lib_install = env.SharedLibrary(join(env['PREFIX'], 'lib', 'ccasm'), 
                                      env['CCASM_SOBJ'], 
                                      CPPPATH=cpppath,
                                      LIBPATH=install_lib_paths, 
                                      LINKFLAGS=linkflags,
                                      LIBS=libs + ['casm'])
Export('ccasm_lib_install')
env.Alias('ccasm_lib_install', ccasm_lib_install)
env['INSTALL_TARGETS'] = env['INSTALL_TARGETS'] + [ccasm_lib_install]

if 'ccasm_lib_install' in COMMAND_LINE_TARGETS:
    env['IS_INSTALL'] = 1


# Include Install instructions
casm_include_install = env.Install(join(env['PREFIX'], 'include'), join(env['INCDIR'], 'casm'))

Export('casm_include_install')
env.Alias('casm_include_install', casm_include_install)
env.Clean('casm_include_install', join(env['PREFIX'], 'include','casm'))
env['INSTALL_TARGETS'] = env['INSTALL_TARGETS'] + casm_include_install

if 'casm_include_install' in COMMAND_LINE_TARGETS:
  env['IS_INSTALL'] = 1

##### Call all SConscript files for executables

# build apps/casm
SConscript(['apps/casm/SConscript'], {'env':env})

# tests/unit
SConscript(['tests/unit/SConscript'], {'env': env})


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
  path = join(env['PREFIX'], 'include', 'casm')
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

##### Configuration checks
if 'configure' in COMMAND_LINE_TARGETS:
  
  def CheckBoost_prefix(conf, boost_prefix):
    conf.Message('BOOST_PREFIX: ' + str(boost_prefix) + '\n')
    conf.Message('Checking for boost headers... ') 
    _path = conf.env.include_path(boost_prefix, 'boost')
    if _path is not None:
      conf.Message('found ' + join(_path, 'boost') + '... ')
      res = 1
    else:
      res = 0
    conf.Result(res)
    return res
  
  def CheckBoost_libname(conf, name):
    if conf.env[name] is None:
      conf.Message("Could not find boost library: " + name)
      res = 0
    else:
      conf.Message("Found boost library: " + i + " as: " + env[i])
      res = 1
    conf.Result(res)
    return res
  
  def CheckBoost_version(conf, version):
    # Boost versions are in format major.minor.subminor
    v_arr = version.split(".")
    version_n = 0
    if len(v_arr) > 0:
        version_n += int(v_arr[0])*100000
    if len(v_arr) > 1:
        version_n += int(v_arr[1])*100
    if len(v_arr) > 2:
        version_n += int(v_arr[2])

    conf.Message('Checking for Boost version >= %s... ' % (version))
    
    ret = conf.TryRun("""
    #include <boost/version.hpp>

    int main() 
    {
        return BOOST_VERSION >= %d ? 0 : 1;
    }
    """ % version_n, '.cpp')[0]
    conf.Result(ret)
    
    if not ret:
      print "Found Boost version:", conf.TryRun("""
      #include <boost/version.hpp>
      #include <iostream>
      int main() 
      {
          std::cout << BOOST_VERSION / 100000 << "."      // major version
                    << BOOST_VERSION / 100 % 1000 << "."  // minor version
                    << BOOST_VERSION % 100                // patch level
                    << std::endl;
          return 0;
      }
      """, ".cpp")[1],
    
    return ret
  
  def CheckBOOST_NO_CXX11_SCOPED_ENUMS(conf):
    conf.Message('Checking if BOOST_NO_CXX11_SCOPED_ENUMS setting is correct... ')
    if ARGUMENTS.get('boost_no_cxx11_scoped_enums', False):
      conf.Message('  Was set via scons command line --boost_no_cxx11_scoped_enums...')
    elif 'CASM_BOOST_NO_CXX11_SCOPED_ENUMS' in os.environ:
      conf.Message('  Was set via environment variable CASM_BOOST_NO_CXX11_SCOPED_ENUMS...')
    
    BOOST_NO_CXX11_SCOPED_ENUMS_test="""
    #include "boost/filesystem.hpp"
    int main(int argc, char* argv[]) {
      boost::filesystem::copy_file("foo", "bar");
      return 0;
    }
    """
    res = conf.TryLink(BOOST_NO_CXX11_SCOPED_ENUMS_test, ".cpp")
    conf.Result(res)
    return res
  
  def check_module(module_name):
    print "Checking for Python module '" + module_name + "'... ",
    try:
      res = imp.find_module(module_name)
      if res:
        print "yes"
        return 1
    except ImportError:
      pass
    for item in sys.path:
      importer = sys.path_importer_cache.get(item)
      if importer:
        try:
          result = importer.find_module(module_name, [item])
          if result:
            print "yes"
            return 1
        except ImportError:
          pass
    print "no"
    return 0
  
  conf = Configure(
    env.Clone(LIBPATH=install_lib_paths,
              CPPPATH=cpppath),
    custom_tests = {
      'CheckBoost_prefix' : CheckBoost_prefix,
      'CheckBoost_libname' : CheckBoost_libname,
      'CheckBoost_version' : CheckBoost_version,
      'CheckBOOST_NO_CXX11_SCOPED_ENUMS': CheckBOOST_NO_CXX11_SCOPED_ENUMS})
  
  def if_failed(msg):
    print "\nConfiguration checks failed."
    print msg
    exit()
  
  # Note: CheckLib with autoadd=1 (default), because some libraries depend on each other
  
  for x in ['z', 'dl']:
    if not conf.CheckLib(env[x]):
      if_failed("Please check your installation")
  if not conf.CheckLib(env['dl']):
    if_failed("Please check your installation")
  if not conf.CheckBoost_prefix(env['BOOST_PREFIX']):
    if_failed("Please check your boost installation or the CASM_BOOST_PREFIX environment variable")
  for i in boost_libs:
    if not conf.CheckBoost_libname(i):
      if_failed("Please check your boost installation or the CASM_BOOST_PREFIX environment variable")
    if not conf.CheckLib(env[i]):
      if_failed("Please check your boost installation or the CASM_BOOST_PREFIX environment variable")
  if not conf.CheckBoost_version('1.54'):
    if_failed("Please check your boost version") 
  if not conf.CheckBOOST_NO_CXX11_SCOPED_ENUMS():
    if_failed("Please check your boost installation or the CASM_BOOST_NO_CXX11_SCOPED_ENUMS environment variable")
  for module_name in ['numpy', 'sklearn', 'deap', 'pandas']:
    if not check_module(module_name):
      if_failed("Python module '" + module_name + "' is not installed")
  if not check_module('pbs'):
      print """
      Python module 'pbs' is not installed
        This module is only necessary for setting up and submitting DFT jobs
        **It is not the module available from pip**"
        It is available from: https://github.com/prisms-center/pbs
      """
  
  
  print "Configuration checks passed."
  exit()
  
