from distutils.core import setup
import glob
setup(name='casm',
      version='supertesting',
      url='git@github.com:goirijo/CASMcode-dev.git',
      description='An interface between CASM and DFT codes.',
      packages=['casm','casm.vaspwrapper', 'casm.learn', 'casm.project', 'casm.plotting'],
      scripts=glob.glob('scripts/*')
      )
