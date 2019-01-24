from __future__ import absolute_import, division, print_function, unicode_literals

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
    return re.sub(' +', ' ', re.sub(chars, '', str(val).strip()))


def confname_as_jobname(configname):
    """
    Returns only the configuration name (SCEL... etc) i.e. for job name in queue
    :param
        configname: Name of the configuration
    :return
        configuration name without path
    """
    return configname.split('/')[-1]
