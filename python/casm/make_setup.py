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
setup(name='casm',
      version='{0}',
      url='{1}',
      description='CASM Python API and interface with DFT codes.',
      packages=['casm','casm.vaspwrapper', 'casm.learn', 'casm.project', 'casm.plotting'],
      scripts=glob.glob('scripts/*')
      )
""".format(version, url)

f = open('setup.py', 'w')
f.write(setup_string)
f.close()
