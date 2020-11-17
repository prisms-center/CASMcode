echo "printenv: "
printenv

echo "$GXX --version: "
$GXX --version

WD=`pwd`
echo "wd: "$WD

# give location of conda-provided g++
ln -s $GCC $PREFIX/bin/gcc
ln -s $GXX $PREFIX/bin/g++
echo "using gcc : : $GXX ;" > tools/build/src/user-config.jam \
  && ./bootstrap.sh \
    --prefix=$PREFIX \
    --with-libraries=system,filesystem,program_options,regex,chrono,timer,test \
  && ./b2 cxxflags="-std=c++17 -O3" -j $CASM_NCPU \
  && ./b2 install -j $CASM_NCPU
