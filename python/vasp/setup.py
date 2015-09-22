from distutils.core import setup
import glob
setup(name='vasp',
      version='bpmaster',
      url='',
      description='A wrapper for submitting and collecting VASP jobs',
      packages=['vasp','vasp.io']
      )
