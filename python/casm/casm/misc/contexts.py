import os
import sys
from mock import patch
from contextlib import contextmanager
from StringIO import StringIO

@contextmanager
def working_dir(wd):
    """
    Usage:
        # do some things in a different working directory, then return to the
        # the current working directory
        with working_dir(wd=some_path):
            ... some code ...
    
    Args:
        wd (str): working directory to use while in context
    """
    orig_wd = os.getcwd()
    if wd is None:
        wd = orig_wd
    os.chdir(wd)
    try:
        yield
    finally:
        os.chdir(orig_wd)

@contextmanager
def captured_output(wd=None):
    """
    Usage:
        # do some things and capture standard output and error
        # as StringIO objects sout, serr
        # optionally, pass 'wd', a different working directory in 
        # which to use while executing the code
        with captured_output(wd=some_path) as (sout, serr):
            ... some code ...
        # print sout and serr nicely:
        print_stringIO(sout)
        print_stringIO(serr)
    
    Args:
        wd (str): working directory to use while in context
    """
    with working_dir(wd):
        new_out, new_err = StringIO(), StringIO()
        old_out, old_err = sys.stdout, sys.stderr
        try:
            sys.stdout, sys.stderr = new_out, new_err
            yield sys.stdout, sys.stderr
        finally:
            sys.stdout, sys.stderr = old_out, old_err

def print_stringIO(strio):
    print "\n----\n", strio.getvalue(), "\n----"

