echo "printenv: "
printenv

echo "g++ --version: "
g++ --version

WD=`pwd`
echo "wd: "$WD

# fix: clang: error: unknown argument: '-fcoalesce-templates'
sed -i.bak -e '141,149 s/^/#/' tools/build/src/tools/darwin.jam

./bootstrap.sh \
    --prefix=$PREFIX \
    --with-libraries=system,filesystem,program_options,regex,chrono,timer,test \
  && ./b2 cxxflags="-std=c++17 -O3" -j $CASM_NCPU \
  && ./b2 install -j $CASM_NCPU
