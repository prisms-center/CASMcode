
CASM C++ library and `ccasm` CLI program
========================================


For Developers
--------------

We recommend developing in the CASM conda development environment. This can be created from the ``CASMcode`` directory by doing:

```
# Location of existing conda installation, or location to newly install conda
export CASM_CONDA_DIR=${_CONDA_ROOT:-"$HOME/.local/conda"}

# Conda environment name
export CASM_ENV_NAME="casm"

# Set location  for CASM test projects to be downloaded from Materials Commons
#  - This is optional, but required for some tests to run
#  - Contact the CASM developers for access and the project ID. You will need
#    to create an account on Materials Commons before we can give you access.
# export CASM_TEST_PROJECTS_DIR="$HOME/.local/CASM_test_projects"
# export CASM_TEST_PROJECTS_ID="<... id here ...>"
# export MC_API_KEY="<... your API key ...>"

bash build_scripts/conda-devel.sh
```

This will download and install conda if it doesn't currently exist, and create a conda environment with the boost libraries required for CASM installed.  On OSX, it will use the system-installed clang compiler, while on linux it will install and use GCC.



Dependencies
------------

### C++11

CASM must be compiled with a compiler that supports the C++11 standard. Testing is done with gcc-7 and Xcode-9.4.

*Note: On Mac OS X it is necessary to have the Xcode command line tools installed. To install, run ``xcode-select --install`` from the Terminal.*

### Boost

The CASM conda development environment includes the subset of the boost library currently used by CASM. Contact the CASM developers if a new feature requires additional boost libraries.

*Important: Boost should be compiled using the same compiler that you will use to compile CASM.*


### Included with CASM

CASM is distributed with the following dependencies:

- **Eigen, v3.1.3**: [http://eigen.tuxfamily.org](http://eigen.tuxfamily.org)

- **JSON Spirit, v4.06**: [http://www.codeproject.com/Articles/20027/JSON-Spirit-A-C-JSON-Parser-Generator-Implemented](http://www.codeproject.com/Articles/20027/JSON-Spirit-A-C-JSON-Parser-Generator-Implemented)

- **Mersenne Twister random number generator -- a C++ class MTRand, v1.1**,  28 September 2009, by Richard J. Wagner, wagnerr@umich.edu, based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus. For more information, see the inventors' web page at [http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html).

- **Qhull, 2015.0.6 2015/10/20**, Barber, C.B., Dobkin, D.P., and Huhdanpaa, H.T., "The Quickhull algorithm for convex hulls," ACM Trans. on Mathematical Software, 22(4):469-483, Dec 1996, http://www.qhull.org.

- **gzstream, v1.5**: [http://www.cs.unc.edu/Research/compgeom/gzstream/](http://www.cs.unc.edu/Research/compgeom/gzstream/)


Build CASM from source
----------------------

With the CASM conda development environment activated, from the "CASMcode" directory, do:

```
bash build.sh
```

Behind the scenes, this will:

1. Run a Python script `make_Makemodule.py` that updates the autoconf and automake files to reflect any new files added to the CASM `include`,  `src`, and `tests` directories (assuming they follow standard naming conventions). See documentation in the script for details.
2. Run autoreconf to generate the `configure` script and execute it.
3. Execute `make`.

On subsequent runs, steps (1) and (2) will be skipped if `configure` exists. To take into account the addition or deletion or files, or to re-configure, do: `rm configure` and re-run.

Options:

- Set `CASM_NCPU=4` to set the `-j` option and compile in parallel
- Set `CASM_CXXFLAGS` for compiler options. The default is:
  - `-O3 -Wall -fPIC --std=c++11 -DNDEBUG -Wno-deprecated-register -Wno-ignored-attributes -Wno-deprecated-declarations`
- See other options in:
  - For OSX: `build_scripts/travis-variables-osx`
  - For Linux: `build_scripts/travis-variables-linux`

To clean up build products:

```
make clean
bash clean.sh
```

Run all tests:
--------------

```
bash build_test.sh
```

Behind the scenes, this will:

1. Do everything that `build.sh` does (build the c++ libraries and programs).
2. Execute `make check` to run c++ tests
3. Run Python tests

Options (beyond those for `build.sh`):

- Set ``CASM_SKIP_CPP_TESTS`` to non-zero length to skip c++ tests.
- Do ``unset CASM_SKIP_CPP_TESTS`` to re-enable c++ tests.
- Set ``CASM_SKIP_PYTHON_TESTS`` to non-zero length to skip Python tests.
- Do ``unset CASM_SKIP_PYTHON_TESTS`` to re-enable Python tests.

See also ``tests/README.md`` for information on writing C++ tests and running select tests.
See also ``python/casm/README.md`` for information on Python tests.

To clean up test products:

```
bash checkclean.sh
```


Install from source:
--------------------

```
bash build_install.sh
```

Uninstall from source:
--------------------

```
bash uninstall.sh
```
