"""test_casm/test_vasp/misc.py"""
import os
from os.path import exists, dirname
import shutil
import distutils.spawn
import warnings

def cp_input(input, output):
    if exists(output):
        shutil.rmtree(output)
    if not exists(dirname(output)):
        os.makedirs(dirname(output))
    shutil.copytree(input, output)

has_vasp = None

def before_all():
    global has_vasp
    has_vasp = (distutils.spawn.find_executable('vasp') is not None)

    if not has_vasp:
        warnings.warn("'vasp' executable not detected: will test behaviour for system without VASP")

    if 'CASM_VASP_POTCAR_DIR' not in os.environ:
        raise Exception("'CASM_VASP_POTCAR_DIR' environment variable not found") 


