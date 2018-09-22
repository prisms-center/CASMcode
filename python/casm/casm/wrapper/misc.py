"""Functions used by wrappers"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import os
import re

def jobname(configname):
    """Return a name for a submitted job for configuration with 'configname'
    
    Args:
    
        configname (str): Name of configuration
    
    Returns: configname.replace(os.sep, '.')
    """
    return configname.replace(os.sep, '.')

def remove_chars(val, chars):
    return re.sub(' +', ' ', re.sub(chars,'',str(val).strip()))
