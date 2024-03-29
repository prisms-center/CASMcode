
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

bash build_scripts/conda-devel.sh
```

This will download and install conda if it doesn't currently exist, and create a conda environment with the boost libraries required for CASM installed.  On OSX, it will use the system-installed clang compiler, while on linux it will install and use GCC.


Using a seperate build directory
--------------------------------

With dependencies installed, CASM can be built using:
```
bash bootstrap.sh
mkdir build && cd build
../configure
make -j4
make install
```

The configuration step may require customization. Here is an example when building in a conda environment:
```
#
CASM_PREFIX=$CONDA_PREFIX

CASM_CXXFLAGS="-O3 -Wall -DNDEBUG -fcolor-diagnostics -I$CASM_PREFIX/include"
CASM_CC="ccache cc"
CASM_CXX="ccache c++"
CASM_PYTHON="python"
CASM_LDFLAGS="-L$CASM_PREFIX/lib -Wl,-rpath,$CASM_PREFIX/lib"
CASM_CONFIGFLAGS="--prefix=$CASM_PREFIX "
CASM_CONFIGFLAGS+="--with-zlib=$CASM_PREFIX "
CASM_CONFIGFLAGS+="--with-boost=$CASM_PREFIX "

../configure CXXFLAGS="${CASM_CXXFLAGS}" CC="$CASM_CC" CXX="$CASM_CXX" PYTHON="$CASM_PYTHON" LDFLAGS="${CASM_LDFLAGS}" ${CASM_CONFIGFLAGS}

```


Dependencies
------------

### C++17

CASM must be compiled with a compiler that supports the C++17 standard. Testing is done with gcc-9 and Xcode-11.6.

*Note: On Mac OS X it is necessary to have the Xcode command line tools installed. To install, run ``xcode-select --install`` from the Terminal.*

### Boost

The CASM conda development environment includes the subset of the boost library currently used by CASM. We are aiming to remove the boost dependency in a future version.

*Important: Boost should be compiled using the same compiler that you will use to compile CASM.*


### Included with CASM

CASM is distributed with the following dependencies:

- **Eigen, v3.1.3**: <http://eigen.tuxfamily.org>

  - Eigen is included with CASM as a git submodule.

- **JSON Spirit, v4.06**: <http://www.codeproject.com/Articles/20027/JSON-Spirit-A-C-JSON-Parser-Generator-Implemented>

- **Mersenne Twister random number generator -- a C++ class MTRand, v1.1**,  28 September 2009, by Richard J. Wagner, wagnerr@umich.edu, based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus. For more information, see the inventors' web page at <http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html>.

- **Qhull, 2015.0.6 2015/10/20**, Barber, C.B., Dobkin, D.P., and Huhdanpaa, H.T., "The Quickhull algorithm for convex hulls," ACM Trans. on Mathematical Software, 22(4):469-483, Dec 1996, http://www.qhull.org.

- **gzstream, v1.5**: <http://www.cs.unc.edu/Research/compgeom/gzstream/>

For testing:

- **googletest**: <https://github.com/google/googletest>

  - googletest is included with CASM as a git submodule


### Not included with CASM ###

- **clang-format**: <https://clang.llvm.org/docs/ClangFormat.html>

  - `clang-format`, with `-style=google`, is used for C++ code formatting. It can be found in most package managers (`apt-get`, `brew`, `yum`, ...). Use the `stylize.sh` script to format files staged for commit. To enforce that it is run before committing, place the `pre-commit` file, or equivalent code, in `.git/hooks`.

- **jq**: <https://stedolan.github.io/jq/>

  - `jq` is used for inspecting and manipulating JSON for testing purposes. It can be found in most package managers.


Build CASM from source
----------------------

Clone CASM, then init and fetch submodules with:

```
git submodule init
git submodule update
```

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
  - `-O3 -Wall -fPIC --std=c++17 -DNDEBUG -Wno-deprecated-register -Wno-ignored-attributes -Wno-deprecated-declarations`
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

The first time only, clone the googletest submodule:

```
cd submodules/googletest
git submodule init
git submodule update
```

Then, you can run:

```
bash build_test.sh
```

Behind the scenes, this will:

1. Do everything that `build.sh` does (build the c++ libraries and programs).
2. Execute `make check` to run c++ tests

Options (beyond those for `build.sh`):

- Set ``CASM_SKIP_CPP_TESTS`` to non-zero length to skip c++ tests.
- Do ``unset CASM_SKIP_CPP_TESTS`` to re-enable c++ tests.

See also ``tests/README.md`` for information on writing C++ tests and running select tests.

To clean up test products:

```
bash checkclean.sh
```

Testing in a seperate build directory
-------------------------------------

If building in a seperate directory, some additional configuration may be necessary before running `make check` for the CASM tests to properly compile generated code. Here is an example when building in a conda environment:
```
export CASM_CXXFLAGS="-O3 -Wall -DNDEBUG -fPIC --std=c++17 -DGZSTREAM_NAMESPACE=gz -I$CONDA_PREFIX/include -I../include"
export CASM_SOFLAGS="-shared -L$CONDA_PREFIX/lib -L.libs -lz"
export CASM_PREFIX=$CONDA_PREFIX
export CASM_BOOST_PREFIX=$CONDA_PREFIX
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
