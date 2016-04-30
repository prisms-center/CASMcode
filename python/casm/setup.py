from distutils.core import setup
import glob
setup(name='casm',
      version='0.2.X_coord_tests',
      url='https://github.com/jcthomas/CASMCode-dev',
      description='An interface between CASM and DFT codes.',
      packages=['casm','casm.vaspwrapper', 'casm.learn', 'casm.project', 'casm.plotting'],
      scripts=glob.glob('scripts/*')
      )
