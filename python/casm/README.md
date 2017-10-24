The CASM Python packages
========================

The `casm` Python packages provide a Python interface to the CASM libraries, implement wrappers to fitting methods and DFT software, and provide other tools for plotting and analysis.

Dependencies
------------

- **Python** Most testing has been done using Python 2.7.5

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


Install from source
-------------------
From ``CASMcode/python/casm`` directory:

	pip install .

This will install:

- required dependencies
- ``casm`` a collection of python packages that provide a Python interface to the CASM libraries, implement wrappers to fitting methods and DFT software, and provide other tools for plotting and analysis.
- ``casm-learn`` a program for fitting effective cluster interactions (ECI)
- ``casm-calc`` a program that helps setup and run high throughput *ab initio* calculations
- ``casm.plot.<type>`` programs to generate ``bokeh`` plots



The functionality provided by ``casm-calc`` is also provided by the legacy scripts:

- ``$CASM_PREFIX/bin/vasp.setup`` a script for setting up VASP jobs
- ``$CASM_PREFIX/bin/vasp.relax`` a script for setting up and submitting VASP jobs
- ``$CASM_PREFIX/bin/vasp.relax.report`` a script for setting up and submitting VASP jobs


Uninstall
---------
 
	pip uninstall casm


Testing dependencies
--------------------

- **pytest** ([https://docs.pytest.org/en/latest/](https://docs.pytest.org/en/latest/))
- **mock** ([https://pypi.python.org/pypi/mock](https://pypi.python.org/pypi/mock))

From ``CASMcode/python/casm`` directory:

	pip install -r test_requirements.txt


Testing configuration
---------------------

Environment variables:

- ``PYTHONPATH``: Set as necessary to find test requirements
- ``CASM_VASP_POTCAR_DIR``: Location of VASP POTCAR files
	- required for some casm.vasp and casm.vaspwrapper tests
- ``CASM_TEST_PROJECTS_DIR``: Location of test CASM projects
	- required for some tests
	- ask the CASM developers for access to the test projects at <casm-developers@lists.engr.ucsb.edu>


Writing Tests
-------------

See documentation for basics on writing tests:

- [pytest documentation](https://docs.pytest.org/en/latest/)
- pytest can run [unittest](https://docs.python.org/2/library/unittest.html) tests, and includes additional features
- See ``test_casm/test_vasp/misc.py`` and ``test_casm/test_vasp/test_relax.py`` for an (incomplete) example implementing and using a base TestCase class that performs common setup for all ``casm.vasp`` tests. This includes checking if the ``vasp`` executable or pseudopotentials exist, and enables using ``skip`` files as described below to provide some additional control over which tests get run.
- Installation of the `casm` CLI and shared libraries is required for testing, but installation of the Python package should not be required for testing


Skip tests
----------

- Add a file named `'skip'` to a test directory to skip all test classes in that directory.
- Add a file named `'skip_TestCase'` to skip the `TestCase` class in that directory.


Run all tests
-------------

From ``CASMcode/python/casm`` directory:

	pytest -r ap -s test_casm


Run all tests in a file
-----------------------

From ``CASMcode/python/casm`` directory:

	pytest -r ap -s test_casm/test_vasp/test_relax.py


Run all tests in a fixture
--------------------------

From ``CASMcode/python/casm`` directory:

	pytest -r ap -s test_casm/test_vasp/test_relax.py::TestRelax


Run a particular test case
--------------------------

From ``CASMcode/python/casm`` directory:

	pytest -r ap -s test_casm/test_vasp/test_relax.py::TestRelax::test_run