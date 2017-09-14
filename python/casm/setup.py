
from distutils.core import setup
import glob
from casm import __version__
setup(name='casm',
      version=__version__,
      url='https://github.com/prisms-center/CASMcode',
      description='CASM Python interface, tools, and wrappers.',
      author='CASM developers',
      author_email='casm-developers@lists.engr.ucsb.edu',
      license='LGPL2.1',
      packages=[
          'casm',
          'casm.misc',
          'casm.vasp',
          'casm.vaspwrapper',
          'casm.seqquest',
          'casm.questwrapper',
          'casm.quantumespresso',
          'casm.qewrapper', 
          'casm.learn',
          'casm.project',
          'casm.plotting'],
      scripts=glob.glob('scripts/*'),
      install_requires=['prisms_jobs', 'scipy', 'pandas', 'scikit-learn', 'bokeh==0.12.3', 'tornado==4.3']
      )
