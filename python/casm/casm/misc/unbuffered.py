"""Write without buffering"""
from __future__ import print_function

def _print(*args, sep=' ', end='\n', file=sys.stdout):
    """Print and always flush"""
    print(*args, sep=sep, end=end, file=file)
    file.flush()
