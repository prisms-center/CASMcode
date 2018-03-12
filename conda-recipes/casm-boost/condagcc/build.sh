echo "printenv: "
printenv

echo "$GXX --version: "
$GXX --version

WD=`pwd`
echo "wd: "$WD

# give location of conda-provided g++
echo "using gcc : : $GXX ;" > tools/build/src/user-config.jam \
  && ./bootstrap.sh \
    --prefix=$PREFIX \
    --with-libraries=system,filesystem,program_options,regex,chrono,timer,test \
  && ./b2 cxxflags="-std=c++11 -O3" -j $NCPUS \
  && ./b2 install -j $NCPUS 

