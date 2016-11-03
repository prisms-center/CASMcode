import sys
if sys.version_info >= (2, 6):
  from setuptools import setup
else:
  from distutils.core import setup
import glob
setup(name='seqquest',
      version='0.01',
      url='https://github.com/prisms-center/CASMcode',
      description='A wrapper for submitting and collecting SeqQuest jobs developed for CASM',
      author='Elizabeth Decolvenaere',
      author_email='edecolv@sandia.gov',
      license='LGPL2.1',
      packages=['seqquest','seqquest.seqquest_io', 'seqquest.seqquest_io.geom', 'seqquest.seqquest_io.lcao_in']
      )
