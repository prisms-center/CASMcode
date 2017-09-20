
CASM C++ library and `casm` CLI program
=======================================

## Dependencies

### Compilation ###

**C++11**

CASM must be compiled with a compiler that supports the C++11 standard. Testing is done with gcc-4.8.5 and clang-800.0.42.1.

**If using Mac OS X - Xcode command-line tools**

On OS X it may be necessary to have the Xcode command line tools installed. To install, run ``xcode-select --install`` from the Terminal.

### SCons

CASM is built using SCons. SCons is available through many package managers, or can be downloaded from the SCons website: [http://www.scons.org](http://www.scons.org).  

**Mac OS X**

Using Macports:

``
sudo port install scons
``

or using Homebrew:

``
brew install scons
``

**Linux Ubuntu**

Install SCons with:

``
sudo apt-get install scons
``

**From the SCons website**

Scons can be downloaded and installed from source following instructions found at the Scons website: [http://www.scons.org/doc/2.3.0/HTML/scons-user/x167.html](http://www.scons.org/doc/2.3.0/HTML/scons-user/x167.html).


### Boost

CASM uses several Boost libraries, which are often available installed on many computing clusters. you can install Boost yourself via a package management tool, or by downloading from the Boost website: [http://www.boost.org](http://www.boost.org). CASM uses the system, filesystem, program_options, and unit_test_framework libraries, and their dependencies. Most CASM testing has been performed with Boost version 1.54 or later. Known bugs in version 1.53 and earlier will cause failures.

*Important: Boost should be compiled using the same compiler that you will use to compile CASM.*

**Mac OS X**
Using Macports:

``
sudo port install boost
``

or using Homebrew:

``
brew install boost
``

**Linux Ubuntu**

Install Boost with:

``
sudo apt-get install libboost-all-dev
``

**From the Boost website**

Boost can be downloaded and installed from source following instructions found at the Boost website: [http://www.boost.org/users/download/](http://www.boost.org/users/download/). 


### Python

CASM includes python modules for automating the submission and analysis of VASP calculations. They have been most extensively tested using Python 2.7.5, and should be compatible with versions 2.x+. (*Note however that for recent versions of SCons, support for Python versions before 2.7 has been deprecated.*) The latest version can be obtained from the Python website: [https://www.python.org](https://www.python.org)

Individual module dependencies include:

- **SciPy** ([https://www.scipy.org](https://www.scipy.org)), which can be obtained using one of the methods described on their website:  [http://www.scipy.org/install.html](http://www.scipy.org/install.html). The particular SciPy packages needed are:
	- **numpy**  ([http://www.numpy.org](http://www.numpy.org))
	- **pandas** ([http://pandas.pydata.org](http://pandas.pydata.org))

- **scikit-learn** ([http://scikit-learn.org](http://scikit-learn.org))

- **deap** ([http://deap.readthedocs.io/en/master/](http://deap.readthedocs.io/en/master/)), the Distributed Evolutionary Algorithm Package, used for genetic algorithms.
	- **scoop** ([http://scoop.readthedocs.io/en/latest/](http://scoop.readthedocs.io/en/latest/)), required for deap. 		

- **pbs** The Python module pbs is used to automate submission and management of PBS batch jobs on a cluster. It can be obtained from its GitHub repository: [https://github.com/prisms-center/pbs](https://github.com/prisms-center/pbs). *Note: This is not the pbs module available for installation via pip.*


### Included with CASM

CASM is distributed with the following dependencies:

- **Eigen, v3.1.3**: [http://eigen.tuxfamily.org](http://eigen.tuxfamily.org) 

- **JSON Spirit, v4.06**: [http://www.codeproject.com/Articles/20027/JSON-Spirit-A-C-JSON-Parser-Generator-Implemented](http://www.codeproject.com/Articles/20027/JSON-Spirit-A-C-JSON-Parser-Generator-Implemented)

- **Mersenne Twister random number generator -- a C++ class MTRand, v1.1**,  28 September 2009, by Richard J. Wagner, wagnerr@umich.edu, based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus. For more information, see the inventors' web page at [http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html).

- **Qhull, 2015.0.6 2015/10/20**, Barber, C.B., Dobkin, D.P., and Huhdanpaa, H.T., "The Quickhull algorithm for convex hulls," ACM Trans. on Mathematical Software, 22(4):469-483, Dec 1996, http://www.qhull.org.

- **gzstream, v1.5**: [http://www.cs.unc.edu/Research/compgeom/gzstream/](http://www.cs.unc.edu/Research/compgeom/gzstream/)


## Installation from source

#### Getting CASM

There are two options for downloading CASM:

1. Fork and clone the repository
    - Use this option if you think you may contribute bug fixes or new features
    - The CASM repository is located at: [https://github.com/prisms-center/CASMcode](https://github.com/prisms-center/CASMcode).  
    - For help see: [https://help.github.com/articles/fork-a-repo/](https://help.github.com/articles/fork-a-repo/).

2. Download as archive
    - Use this option if you do not plan on contributing bug fixes or new features
    - CASM can be downloaded as a .zip or .tar.gz archive from: [https://github.com/prisms-center/CASMcode/releases](https://github.com/prisms-center/CASMcode/releases)

#### Configuration
CASM is built using SCons, but some configuration of environment variables may be necessary first. Help can be obtained via:

    scons -h

and is also reproduced here:

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
      
      $CASM_PREFIX
      $CASM_INCLUDEDIR (=$CASM_PREFIX/include)
      $CASM_LIBDIR (=$CASM_PREFIX/lib)
      $CASM_BINDIR (=$CASM_PREFIX/bin)
      $CASM_PYTHON_PREFIX (=$CASM_PREFIX)
        Where to install CASM. By default, this uses 'CASM_PREFIX=/usr/local'. Header 
        files are installed in '$CASM_INCLUDEDIR', shared libraries in '$CASM_LIBDIR', 
        executables in '$CASM_BINDIR', and $CASM_PYTHON_PREFIX is used for the 
        setup.py --prefix option for installing python packages.
      
      $CASM_BOOST_PREFIX
      $CASM_BOOST_INCLUDEDIR (=$CASM_BOOST_PREFIX/include)
      $CASM_BOOST_LIBDIR (=$CASM_BOOST_PREFIX/lib)
        Search path for Boost. '$CASM_BOOST_INCLUDEDIR' is searched for header files, and
        '$CASM_BOOST_LIBDIR' for libraries. 
        Boost and CASM should be compiled with the same compiler.

      $CASM_OPTIMIZATIONLEVEL:
        Sets the -O optimization compiler option. If not set, uses -O3.

      $CASM_DEBUGSTATE:
        Sets to compile with debugging symbols. In this case, the optimization level gets 
        set to -O0, and NDEBUG does not get set.

      $CASM_BOOST_NO_CXX11_SCOPED_ENUMS:
        If defined, will compile with -DCASM_BOOST_NO_CXX11_SCOPED_ENUMS. Use this
        if linking to boost libraries compiled without c++11.
      
      $CASM_BASH_COMPLETION_DIR:
        If defined, bash-completion scripts for CASM will be installed in the 
        location given. If not defined, standard locations will be searched for
        'bash_completion'. 
      
      
      Additional options that override environment variables:
      
      Use 'cxx=X' to set the C++ compiler. Default is chosen by scons.
          'opt=X' to set optimization level, '-OX'. Default is 3.
          'debug=X' with X=0 to use '-DNDEBUG', 
             or with X=1 to set debug mode compiler options '-O0 -g -save-temps'.
             Overrides $CASM_DEBUGSTATE.
          'casm_prefix=X' to set installation directory. Default is '/usr/local'. Overrides $CASM_PREFIX.
          'boost_prefix=X' set boost search path. Overrides $CASM_BOOST_PPREFIX.
          'boost_no_cxx11_scoped_enums=1' to use '-DBOOST_NO_CXX11_SCOPED_ENUMS'.
             Overrides $CASM_BOOST_NO_CXX11_SCOPED_ENUMS.
     
      Use scons -H for help about command-line options.


The script ``casmenv.sh`` provides a list of environment variables that you are recogized by CASM during installation and use.  A copy of this file can be used to configure your environment before installing or using CASM. For instance:
	
	mkdir $HOME/modules
	cp /path/to/CASMcode/casmenv.sh $HOME/modules/casm
	... edit $HOME/modules/casm ...

Then to set your environment before installing or using CASM:
	
	source $HOME/modules/casm

If you are working in a shared computing environment where other modules such as VASP must be imported, they can also be imported in this script.

After setting up your environment, run:

	cd /path/to/CASMcode
	scons configure

to perform a number of configuration checks. Once they pass, you are ready to install CASM.

#### Build and install


Once any necessary environment variables are set, you are ready to build and install. Move to the directory in which the CASM source resides and run ``scons install``:

    cd /path/to/CASMcode
    scons install

*Note: Use 'scons install -j N', where N is number of jobs, to enable multi-threaded compilation. This can make a nice difference.*

This will compile and install:

- ``$CASM_PREFIX/bin/casm`` the primary CASM program
- ``$CASM_PREFIX/bin/casm-learn`` a program for fitting effective cluster interactions (ECI)
- ``$CASM_PREFIX/bin/casm-calc`` a program that helps setup and run high throughput *ab initio* calculations
- ``$CASM_PREFIX/include/casm/`` headers files for ``libcasm``
- ``$CASM_PREFIX/lib/libcasm.*`` a shared library containing much of the implementation. May be ``libcasm.dylib`` on Mac OS X.
- ``$CASM_PREFIX/lib/libccasm.*`` a shared library providing a C interface to ``libcasm.*`` used by the ``casm`` Python package
- ``$CASM_PREFIX/lib/pythonX.Y/site-packages/casm`` a python package that provides an interface between ``casm`` and the software used to calculate training data. Currently only VASP is supported.
- ``$CASM_PREFIX/lib/pythonX.Y/site-packages/vasp`` a python package for running VASP calculations.
- (optional) ``$CASM_PREFIX/bin/casm-complete`` implements bash-completion

The functionality provided by ``casm-calc`` is also provided by the legacy scripts:

- ``$CASM_PREFIX/bin/vasp.setup`` a script for setting up VASP jobs
- ``$CASM_PREFIX/bin/vasp.relax`` a script for setting up and submitting VASP jobs
- ``$CASM_PREFIX/bin/vasp.relax.report`` a script for setting up and submitting VASP jobs




#### Checking installation ####

If ``casm`` is installed correctly, execute ``casm`` from any directory and you will see the ``casm`` help menu:

	*** casm usage ***
	
	casm [--version] <subcommand> [options] [args]
	
	available subcommands:
	  bset
	  composition
	  enum
	  fit
	  format
	  import
	  init
	  perturb
	  query
	  ref
	  run
	  select
	  settings
	  status
	  super
	  sym
	  update
	
	For help using a subcommand use: 'casm <subcommand> --help'
	
	For step by step help use: 'casm status -n'



**Frequently encountered issues**:

- **I tried to install (``scons install``) or uninstall (``scons install -c``), but get errors about not having permission.**
  - If you don't have permissions to write to ``/usr/local/``, as is usual on a computer cluster, you can change the environment variable ``$CASM_PREFIX`` in your configure script ``$HOME/modules/casm`` to a location that you do have permission to use. For instance, you can create a software directory in your home directory:
  
      	cd ~
      	mkdir software
    
    To make the changes take effect open a new session and 
    
    	source $HOME/modules/casm
    
    Then try installing again 
    	
    	cd /path/to/CASMcode
    	scons install

  - If you have administrative access you can install using ``sudo``, although this is not recommended. For example: ``sudo scons install``.

        

- **I installed CASM without errors, but when I try to execute ``casm`` I get the error**:

        $ casm
		-bash: casm: command not found
		
  If ``scons install`` ran without error, this means that ``casm`` was installed in a directory that is not in your $PATH. You can check what directories are searched for executables using ``echo $PATH``. One solution is to:
  1. Uninstall CASM and re-install it in a directory that is in your ``$PATH``. This can be accomplished by uninstalling (``scons install -c``) and changing ``$CASM_PREFIX`` in ``$HOME/modules/casm`` such that ``$CASM_PREFIX/bin`` is in your ``$PATH``. To make the changes take effect open a new session and 
    
    		source $HOME/modules/casm 
   

- **I installed CASM without errors, but when I try to execute ``casm`` I get the error**: 
    
        $ casm
        casm: error while loading shared libraries: libcasm.so: cannot open 
        shared object file: No such file or directory
  
  This means ``casm`` has been installed correctly but the shared library ``libcasm.so`` (or ``libcasm.dylib`` on Mac OS X) is not being found in any of the places being searched.  Possible solutions include:
  - Check that the ``LD_LIBRARY_PATH`` environment variable in ``$HOME/modules/casm`` is specifyign which directory to search for ``libcasm``:
       
        export LD_LIBRARY_PATH=$CASM_PREFIX/lib:$LD_LIBRARY_PATH
        
    To make the changes take effect open a new session and 
    
    	source $HOME/modules/casm
    
    *Note: On Mac OS X, use ``DYLD_FALLBACK_LIBRARY_PATH`` instead of ``LD_LIBRARY_PATH``.*
  - (Linux): Update the default library search path using ``ldconfig``. For example, see [this](http://codeyarns.com/2014/01/14/how-to-add-library-directory-to-ldconfig-cache/).
 
#### For developers: Testing new features ####

If you are developing new features you can run all unit and integration tests from the main repository directory via:

    scons unit

Individual tests can also be run via:

    scons unit
    scons casm_test
    scons eci_search_test
    # replace UnitTestName with a particular unit test:
    # (UnitTestName in tests/unit/*/UnitTestName_test.cpp)
    scons UnitTestName

Individual tests may be cleaned by re-running with any of the above commands with an added ``-c``. For instance ``scons Clexulator -c`` or ``scons casm_test -c``. 

New unit tests using the Boost unit test framework can be added and run by placing a ``UnitTestName_test.cpp`` file in any subdirectory of ``tests/unit``, following the template of existing unit tests. If the unit test needs to link to libraries or creates any files that it doesn't remove by itself, the ``tests/unit/SConscript`` file should be edited to enable link needed libraries and clean generated them.


