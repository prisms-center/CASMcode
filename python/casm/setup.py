from distutils.core import setup
import glob
setup(name='casm',
      version='0.1.0',
      url='https://github.com/prisms-center/CASMcode.git',
      description='An interface between CASM and DFT codes.',
      packages=['casm','casm.vaspwrapper'],
      scripts=glob.glob('scripts/*')
      )
