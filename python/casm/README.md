The CASM Python packages
========================

The `casm-python` Python packages provide a Python interface to the CASM libraries, implement wrappers to fitting methods and DFT software, and provide other tools for plotting and analysis.


Python testing dependencies
--------------------

- **pytest** ([https://docs.pytest.org/en/latest/](https://docs.pytest.org/en/latest/))
- **mock** ([https://pypi.python.org/pypi/mock](https://pypi.python.org/pypi/mock))

From ``CASMcode/python/casm`` directory:

	pip install -r test_requirements.txt


Testing configuration
---------------------

Environment variables control whether some tests run:

- ``CASM_VASP_POTCAR_DIR``: Location of VASP POTCAR files
	- required for some `casm.vasp` and `casm.vaspwrapper` tests
- ``CASM_TEST_PROJECTS_DIR``: Location of test CASM projects
	- required for some tests
	- ask the CASM developers for access to the test projects at <casm-developers@lists.engr.ucsb.edu>


Writing tests
-------------

See documentation for basics on writing tests:

- [pytest documentation](https://docs.pytest.org/en/latest/)
- pytest can run [unittest](https://docs.python.org/2/library/unittest.html) tests, and includes additional features
- See ``test_casm/test_vasp/misc.py`` and ``test_casm/test_vasp/test_relax.py`` for an (incomplete) example implementing and using a base TestCase class that performs common setup for all ``casm.vasp`` tests. This includes checking if the ``vasp`` executable or pseudopotentials exist, and enables using ``skip`` files as described below to provide some additional control over which tests get run.
- Installation of the `casm` CLI and shared libraries is required for testing, but installation of the Python package should not be required for testing


Running tests
-------------

From ``CASMcode`` directory:

Run all tests:

```
bash build_test.sh
```

Behind the scenes, this will:

1. Do everything that `build.sh` does (build the C++ libraries and programs).
2. Execute `make check` to run C++ tests
3. Run Python tests

Options:

- Set ``CASM_SKIP_CPP_TESTS`` to non-zero length to skip C++ tests.
- Do ``unset CASM_SKIP_CPP_TESTS`` to re-enable C++ tests.
- Set ``CASM_SKIP_PYTHON_TESTS`` to non-zero length to skip Python tests.
- Do ``unset CASM_SKIP_PYTHON_TESTS`` to re-enable Python tests.
- Use the ``CASM_PYTEST_ARGS`` environment variable to pass arguments to pytest. Default is `test_casm`.
- Use the ``CASM_PYTEST_OPTIONS`` environment variable to set pytest options. Default is ``"-r ap -s"``.
- Add a file named `'skip'` to a test directory to skip all test classes in that directory.
- Add a file named `'skip_TestCase'` to skip the `TestCase` class in that directory.


Examples running tests:
-----------------------

Run all tests in a file:

```
export CASM_PYTEST_ARGS="test_casm/test_vasp/test_relax.py"
bash build_test.sh
```

Run all tests in a fixture:

```
export CASM_PYTEST_ARGS="test_casm/test_vasp/test_relax.py::TestRelax"
bash build_test.sh
```

Run a particular test case:

```
export CASM_PYTEST_ARGS="test_casm/test_vasp/test_relax.py::TestRelax::test_run"
bash build_test.sh
```


Dependencies
------------

Dependencies should be installed automatically. The primary dependencies are:

- **Python** Currently, we are testing using Python 3.6.  We are maintaining Python 2.7 compatibility for the near future.

Individual Package dependencies include:

- **SciPy** ([https://www.scipy.org](https://www.scipy.org)), which can be obtained using one of the methods described on their website:  [http://www.scipy.org/install.html](http://www.scipy.org/install.html). The particular SciPy packages needed are:
	- **numpy**  ([http://www.numpy.org](http://www.numpy.org))
	- **pandas** ([http://pandas.pydata.org](http://pandas.pydata.org))
	- **bokeh** ([https://bokeh.pydata.org/en/latest/](https://bokeh.pydata.org/en/latest/))

- **scikit-learn** ([http://scikit-learn.org](http://scikit-learn.org))

- **deap** ([http://deap.readthedocs.io/en/master/](http://deap.readthedocs.io/en/master/)), the Distributed Evolutionary Algorithm Package, used for genetic algorithm

- **prisms_jobs** ([https://prisms-center.github.io/prisms_jobs_docs](https://prisms-center.github.io/prisms_jobs_docs))


Generating html documentation
-----------------------------
From ``CASMcode/python/casm`` directory:

	# Install sphinx requirements
	pip install -r doc_requirements.txt

	# Generate an index.rst file including all casm subpackages
	python build_doc_api_index.py

	# Generate docs
	python setup.py build_sphinx

	# Open
	open doc/build/html/index.html
