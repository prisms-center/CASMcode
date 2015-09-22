from distutils.core import setup
import glob
setup(name='casm',
      version='bpmaster',
      url='',
      description='An interface between CASM and DFT codes.',
      packages=['casm','casm.vaspwrapper'],
      scripts=glob.glob('scripts/*')
      )
