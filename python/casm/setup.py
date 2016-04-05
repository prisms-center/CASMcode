from distutils.core import setup
import glob
setup(name='casm',
      version='0.1.1_strain_enum',
      url='https://github.com/jcthomas/CASM_JCT.git',
      description='An interface between CASM and DFT codes.',
      packages=['casm','casm.vaspwrapper'],
      scripts=glob.glob('scripts/*')
      )
