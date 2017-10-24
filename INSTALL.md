
CASM C++ library and `casm` CLI program
=======================================

*Note: This only gives installation instructions for the `casm` command line program and c/c++ libraries. To install the CASM Python packages and programs see ``python/casm/README.md``.*
 

Dependencies
------------

### C++11 

CASM must be compiled with a compiler that supports the C++11 standard. Testing is done with gcc-4.8.5 and Apple LLVM version 7.3.0 (clang-703.0.31).

*Note: On Mac OS X it may be necessary to have the Xcode command line tools installed. To install, run ``xcode-select --install`` from the Terminal.*

### Boost

CASM uses several Boost libraries, which are often available installed on many computing clusters. you can install Boost yourself via a package management tool, or by downloading from the Boost website: [http://www.boost.org](http://www.boost.org). Most CASM testing has been performed with Boost version 1.54 or later. Known bugs in version 1.53 and earlier will cause failures.

*Important: Boost should be compiled using the same compiler that you will use to compile CASM.*

**Mac OS X**
Using Homebrew:

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


#### Configure, Build, and install

If installing from an archive distribution a ``configure`` script will already be present. If installing from the git repository, first generate the ``configure`` script using:

    ./bootstrap.sh

On many systems, all that is required is:

    ./configure && make && make install

Some systems may require special options for configuration. Help can be obtained from ``./configure --help``. In particular, you many need to use:

- ``CXXFLAGS="-O3 -DNDEBUG -DBOOST_NO_CXX11_SCOPED_ENUMS -Wno-deprecated-register -Wno-deprecated-declarations"``
  - ``-O3`` for maximum optimization
  - ``-DNDEBUG`` to disable debugging mode
  - ``-DBOOST_NO_CXX11_SCOPED_ENUMS`` if boost was compiled without c++11
  - ``-Wno-deprecated-register -Wno-deprecated-declarations`` to disable some compiler warnings
- ``--prefix=PREFIX`` to set the installation location
- ``--with-boost=PATH`` to find boost in a non-standard location
- ``--with-boost-libdir=PATH`` to find boost libraries in a non-standard location
- ``--with-bash-completion-dir=PATH`` to specify where to install the CASM ``bash-completion`` script



This will compile and install:

- ``casm`` the primary CASM program
- CASM headers files
- ``libcasm`` a shared library containing much of the implementation
- ``libccasm`` a shared library providing a C interface to ``libcasm`` used by the ``casm`` Python package
- (optional) ``casm-complete`` a program that implements bash-completion and an associated bash-completion script

#### Checking installation ####

If ``casm`` is installed correctly, execute ``casm`` from any directory and you will see the ``casm`` help menu:

```
-- casm usage -- 

casm [--version] <command> [options] [args]

available commands:
  bset
  composition
  enum
  ...
```
 
#### For developers: Testing new features ####

See ``tests/README.md``
