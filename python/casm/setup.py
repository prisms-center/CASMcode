from setuptools import setup, find_packages
import glob
import os
from casm import __version__

# get console_scripts
def script_str(file):
    name = os.path.splitext(os.path.split(file)[1])[0]
    if name in ['casm_calc', 'casm_learn', 'casm_plot']:
        return name.replace('_','-') + '=casm.scripts.' + name + ':main'
    else:
        return name.replace('_','.') + '=casm.scripts.' + name + ':main'
console_scripts = [script_str(x) for x in glob.glob('casm/scripts/*.py') if x != 'casm/scripts/__init__.py']
print(console_scripts)

setup(name='casm-python',
      version=__version__,
      url='https://github.com/prisms-center/CASMcode',
      description='CASM Python interface, tools, and wrappers.',
      author='CASM developers',
      author_email='casm-developers@lists.engr.ucsb.edu',
      license='LGPL2.1+',
      packages=find_packages(),
      entry_points={
          'console_scripts': console_scripts
      },
      install_requires=[
          'bokeh',
          'deap',
          'mock', 
          'pandas',
          'prisms-jobs',
          'scikit-learn',
          'scipy',
          'sh',
          'tornado'
      ],
      classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering'
      ],
      data_files = [('', ['LICENSE'])]
      )
