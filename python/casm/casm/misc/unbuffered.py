"""Write without buffering"""
from __future__ import print_function

def _print(*args, **kwargs):
    """Print and always flush"""
    sep = kwargs.get('sep', ' ')
    end = kwargs.get('end', '\n')
    file = kwargs.get('file', sys.stdout)
    print(*args, sep=sep, end=end, file=file)
    file.flush()
