echo "printenv: "
printenv

echo "$GXX --version: "
$GXX --version

echo "$CC --version: "
$CC --version

echo "$CXX --version: "
$CXX --version

WD=`pwd`
echo "wd: "$WD

CASM_BASH_COMPLETION_DIR="$PREFIX"/.bash_completion.d
mkdir -p $CASM_BASH_COMPLETION_DIR \
  && printf "for bcfile in $CASM_BASH_COMPLETION_DIR/* ; do\n  . \$bcfile\ndone" > $PREFIX/.bash_completion

cd $WD \
  && pip install --no-cache-dir six \
  && python make_Makemodule.py \
  && ./bootstrap.sh \
  && ./configure \
    CXXFLAGS="-O3 -DNDEBUG -Wno-deprecated-register -Wno-ignored-attributes -Wno-deprecated-declarations" \
    --prefix=$PREFIX \
    --with-zlib=$PREFIX \
    --with-boost-libdir=$PREFIX/lib \
    --with-bash-completion-dir=$CASM_BASH_COMPLETION_DIR \
  && make -j $NCPUS \
  && make install
