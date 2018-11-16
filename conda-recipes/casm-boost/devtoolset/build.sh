echo "printenv: "
printenv

echo "g++ --version: "
g++ --version

WD=`pwd`
echo "wd: "$WD

./bootstrap.sh \
    --prefix=$PREFIX \
    --with-libraries=system,filesystem,program_options,regex,chrono,timer,test \
  && ./b2 cxxflags="-std=c++11 -O3" -j $CASM_NCPU \
  && ./b2 install -j $CASM_NCPU
