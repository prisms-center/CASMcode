import sys

if len(sys.argv) != 3:
  print "Error running casm/make_setup.py."

version = sys.argv[1]
url = sys.argv[2]

setup_string = \
"""import sys
if sys.version_info >= (2, 6):
  from setuptools import setup
else:
  from distutils.core import setup
import glob
setup(name='vasp',
      version='{0}',
      url='{1}',
      description='A wrapper for submitting and collecting VASP jobs developed for CASM',
      packages=['vasp','vasp.io']
      )
""".format(version, url)

f = open('setup.py', 'w')
f.write(setup_string)
f.close()
