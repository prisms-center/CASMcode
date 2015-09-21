
## Dependencies

### Compilation ###

**C++11**

CASM must be compiled with a compiler that supports the C++11 standard.

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

CASM uses several Boost libraries, which are often available installed on many computing clusters. you can install Boost yourself via a package management tool, or by downloading from the Boost website: [http://www.boost.org](http://www.boost.org). CASM uses the system, filesystem, program_options, and unit_test_framework libraries, and their dependencies. Most CASM testing has been performed with Boost version 1.54 or later.

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

**NumPy**

Individual module dependencies include NumPy ([http://www.numpy.org](http://www.numpy.org)), which can be obtained by installing SciPy using one of the methods described on their website:  [http://www.scipy.org/install.html](http://www.scipy.org/install.html).

**pbs**

The Python module pbs is used to automate submission and management of PBS batch jobs on a cluster. It can be obtained from its GitHub repository: [https://github.com/prisms-center/pbs](https://github.com/prisms-center/pbs).


### Included with CASM

CASM is distributed with the following dependencies:

- **Eigen, v3.1.3**: [http://eigen.tuxfamily.org](http://eigen.tuxfamily.org) 

- **JSON Spirit, v4.06**: [http://www.codeproject.com/Articles/20027/JSON-Spirit-A-C-JSON-Parser-Generator-Implemented](http://www.codeproject.com/Articles/20027/JSON-Spirit-A-C-JSON-Parser-Generator-Implemented)

- **Mersenne Twister random number generator -- a C++ class MTRand, v1.1**,  28 September 2009, by Richard J. Wagner, wagnerr@umich.edu, based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus. For more information, see the inventors' web page at [http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html).


## Installation

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
        Search path for dynamic libraries, may need $CASMBOOST_PATH/lib 
        and $CASMPREFIX/lib added to it.
        On Mac OS X, this variable is $DYLD_FALLBACK_LIBRARY_PATH.
        This should be added to your ~/.bash_profile (Linux) or ~/.profile (Mac).
      
      
      Additional options that override environment variables:
      
      Use 'cxx=X' to set the C++ compiler. Default is chosen by scons.
          'opt=X' to set optimization level, '-OX'. Default is 3.
          'debug=X' with X=0 to use '-DNDEBUG', 
                    or with X=1 to set debug mode compiler options '-O0 -g -save-temps'.
                    Overrides $DEBUGSTATE.
          'prefix=X' to set installation directory. Default is '/usr/local'. 
                    Overrides $CASMPREFIX.
          'boost_path=X' set boost search path. Overrides $CASMBOOST_PATH.



For example, on a cluster where Boost is installed in a shared directory ``/home/software/rhel6/boost/1.54.0-gcc-4.7.0`` (*Important: Boost and CASM should be compiled with the same compiler.*), and your executables and Python modules are stored in your userspace at ``$HOME/software``, you could add the following to the ``.bash_profile`` file in your home directory:

    export CASMBOOST_PATH=/home/software/rhel6/boost/1.54.0-gcc-4.7.0
    export CASMPREFIX=$HOME/software
    export PATH=$PATH:$CASMPREFIX/bin
    export CPATH=$CPATH:$CASMPREFIX/include
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CASMPREFIX/lib

and then run ``source $HOME/.bash_profile`` for these changes to take effect. These commands to set environment variables could also be run from the command line, in which case the environment variables only persist until your terminal session is closed. *Note: On Mac OS X, use ``DYLD_FALLBACK_LIBRARY_PATH`` instead of ``LD_LIBRARY_PATH``.*

#### Build and install


Once any necessary environment variables are set, you are ready to build and install. Move to the directory in which the CASM source resides and run ``scons install``:

    cd /path/to/CASM
    scons install

*Note: Use 'scons install -j N', where N is number of jobs, to enable multi-threaded compilation. This can make a nice difference.*

This will compile and install:

- ``$CASMPREFIX/bin/casm`` the primary CASM program
- ``$CASMPREFIX/bin/eci_search`` a program for fitting effective cluster interactions (ECI)
- ``$CASMPREFIX/bin/vasp.setup`` a script for setting up VASP jobs
- ``$CASMPREFIX/bin/vasp.relax`` a script for setting up and submitting VASP jobs
- ``$CASMPREFIX/bin/kpoint_converge`` a script for performing k-point convergence
- ``$CASMPREFIX/include/casm/`` headers files for ``libcasm``
- ``$CASMPREFIX/lib/libcasm.so`` a shared library containing much of the implementation. May be ``libcasm.dylib`` on Mac OS X.
- ``$CASMPREFIX/lib/pythonX.Y/site-packages/casm`` a python package that provides an interface between ``casm`` and the software used to calculate training data. Currently only VASP is supported.
- ``$CASMPREFIX/lib/pythonX.Y/site-packages/vasp`` a python package for running VASP calculations.


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

- I tried to install (``scons install``) or uninstall (``scons install -c``), but get errors about not having permission.
  - If you don't have permissions to write to ``/usr/local/``, as is usual on a computer cluster, you can change the environment variable ``$CASMPREFIX`` to a location that you do have permission to use. For instance, you can create a software directory in your home directory:
  
      	cd ~
      	mkdir software
    
    Then you can edit the ``.bash_profile`` file in your home directory to set your ``$PATH`` and libary search path to include your software directory by adding the lines:
    
        export CASMPREFIX=$HOME/software
        export PATH=$PATH:$CASMPREFIX/bin
        export CPATH=$CPATH:$CASMPREFIX/include        
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CASMPREFIX/lib
        
    To make the changes take effect use ``source ~/.bash_profile`` or open a new session. Then try installing again (``cd /path/to/CASM; scons install``). *Note: On Mac OS X, use ``DYLD_FALLBACK_LIBRARY_PATH`` instead of ``LD_LIBRARY_PATH``.*

  - If you have administrative access you can install using ``sudo``, although this is not recommended. For example: ``sudo scons install``.

        

- I installed CASM without errors, but when I try to execute ``casm`` I get the error:

        $ casm
		-bash: casm: command not found
		
  If ``scons install`` ran without error, this means that ``casm`` was installed in a directory that is not in your $PATH. You can check what directories are searched for executables using ``echo $PATH``. Possible solutions include:
  1. Uninstall CASM and re-install it in a directory that is in your ``$PATH``. This can be accomplished by uninstalling (``scons install -c``) and changing ``$CASMPREFIX`` such that ``$CASMPREFIX/bin`` is in your ``$PATH`` (``export $CASMPREFIX=/some/place/in_my_path/)
  1. Append the location where ``casm`` was installed to your ``$PATH``. This can be accomplished by adding the line ``export PATH=$PATH:/path/to/bin`` to the file ``.bash_profile`` in your home directory, where ``/path/to/bin`` is replaced with the actual path to the directory where ``casm`` is installed. 
   

- I installed CASM without errors, but when I try to execute ``casm`` I get the error: 
    
        $ casm
        casm: error while loading shared libraries: libcasm.so: cannot open 
        shared object file: No such file or directory
  
  This means ``casm`` has been installed correctly but the shared library ``libcasm.so`` (or ``libcasm.dylib`` on Mac OS X) is not being found in any of the places being searched.  Possible solutions include:
  - (Linux): Update the default library search path using ``ldconfig``. For example, see [this](http://codeyarns.com/2014/01/14/how-to-add-library-directory-to-ldconfig-cache/).
  - Change the ``LD_LIBRARY_PATH`` environment variable to specify which directory to search for ``libcasm`` by editing the ``.bash_profile`` file in your home directory to include:
       
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH
        
  *Note: On Mac OS X, use ``DYLD_FALLBACK_LIBRARY_PATH`` instead of ``LD_LIBRARY_PATH``.*

#### For developers: Testing new features ####

If you are developing new features you can run all unit and integration tests from the main repository directory via:

    scons test

Individual tests can also be run via:

    scons unit
    scons casm_test
    scons eci_search_test
    # replace UnitTestName with a particular unit test (UnitTestName in tests/unit/*/UnitTestName_test.cpp)
    scons UnitTestName

Individual tests may be cleaned by re-running with any of the above commands with an added ``-c``. For instance ``scons Clexulator -c`` or ``scons casm_test -c``. In particular, ``scons test`` and ``scons casm_test`` must be cleaned before re-running or there will be errors about trying to initialize a CASM project in a directory where one already exists. 

New unit tests using the Boost unit test framework can be added and automatically run by placing a ``UnitTestName_test.cpp`` file in any subdirectory of ``tests/unit`` by following the template of existing unit tests. If the unit test creates any files that it doesn't remove by itself, the ``tests/unit/SConscript`` file should be edited to enable cleaning them.


