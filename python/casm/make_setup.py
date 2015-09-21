import subprocess

# use git branch for version
version = subprocess.Popen('git rev-parse --abbrev-ref HEAD'.split(), stdout=subprocess.PIPE).communicate()[0].rstrip()

if version[0] == 'v':
  version = version[1:]

# get remote url
url = subprocess.Popen('git config --get remote.origin.url'.split(), stdout=subprocess.PIPE).communicate()[0].rstrip()

setup_string = \
"""from distutils.core import setup
import glob
setup(name='casm',
      version='{0}',
      url='{1}',
      description='An interface between CASM and DFT codes.',
      packages=['casm','casm.vaspwrapper'],
      scripts=glob.glob('scripts/*')
      )
""".format(version, url)

f = open('setup.py', 'w')
f.write(setup_string)
f.close()
