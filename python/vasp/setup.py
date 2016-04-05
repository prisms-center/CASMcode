from distutils.core import setup
import glob
setup(name='vasp',
      version='0.1.1_strain_enum',
      url='https://github.com/jcthomas/CASM_JCT.git',
      description='A wrapper for submitting and collecting VASP jobs',
      packages=['vasp','vasp.io']
      )
