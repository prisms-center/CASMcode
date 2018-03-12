echo "printenv: "
printenv

echo "g++ --version: "
g++ --version

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
    --with-boost=$PREFIX \
    --with-bash-completion-dir=$CASM_BASH_COMPLETION_DIR \
  && make -j $NCPUS \
  && make install
