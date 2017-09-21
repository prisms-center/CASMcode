g++ --version

export CASM_PREFIX=$PREFIX
export CASM_BOOST_PREFIX=$PREFIX
export CASM_BOOST_NO_CXX11_SCOPED_ENUMS=1
scons configure
scons -j 4
scons install -j 4

cd python/casm
pip install .
  
