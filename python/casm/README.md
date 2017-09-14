Testing configuration
---------------------

Environment variables:

- CASM_VASP_POTCAR_DIR: Location of VASP POTCAR files
	- required for casm.vasp, casm.vaspwrapper 


Run all tests
-------------

From ``CASMcode/python/casm`` directory:

	pytest -r ap -s test_casm

Run all tests in a file
-----------------------

From ``CASMcode/python/casm`` directory:

	pytest -r ap -s test_casm/test_vasp/test_relax.py

Run all tests in a class
------------------------

From ``CASMcode/python/casm`` directory:

	pytest -r ap -s test_casm/test_vasp/test_relax.py::TestRelax

Run a particular test
---------------------

From ``CASMcode/python/casm`` directory:

	pytest -r ap -s test_casm/test_vasp/test_relax.py::TestRelax::test_run