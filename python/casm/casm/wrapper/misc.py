"""Functions used by wrappers"""
import os

def jobname(configname):
    """Return a name for a submitted job for configuration with 'configname'
    
    Args:
    
        configname (str): Name of configuration
    
    Returns: configname.replace(os.sep, '.')
    """
    return configname.replace(os.sep, '.')
