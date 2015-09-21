from distutils.core import setup
import glob
setup(name='vasp',
      version='0.1.0',
      url='https://github.com/prisms-center/CASMcode.git',
      description='A wrapper for submitting and collecting VASP jobs',
      packages=['vasp','vasp.io']
      )
