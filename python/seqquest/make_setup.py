import sys

if len(sys.argv) != 3:
  print "Error running casm/make_setup.py."

version = sys.argv[1]
url = sys.argv[2]

setup_string = \
"""
from distutils.core import setup
import glob
setup(name='seqquest',
      version='{0}',
      url='{1}',
      description='A wrapper for submitting and collecting SeqQuest jobs developed for CASM',
      author='CASM developers',
      author_email='casm-developers@lists.engr.ucsb.edu',
      license='LGPL2.1',
      packages=['seqquest','seqquest.seqquest_io']
      )
""".format(version, url)

f = open('setup.py', 'w')
f.write(setup_string)
f.close()
