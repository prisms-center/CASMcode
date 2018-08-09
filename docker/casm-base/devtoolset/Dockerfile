FROM centos:7

LABEL description="CASM build & run container, using devtoolset"

# Make it easier to customize build by allowing build arguments
ARG PYTHON_VERSION
ARG DEVTOOLSET_VERSION

ENV CASM_PYTHON_VERSION=$PYTHON_VERSION
ENV CASM_CONDA_DIR=/home/casmuser/.local/conda
ENV CASM_DEVTOOLSET=devtoolset-$DEVTOOLSET_VERSION


### Stage 1: Install build tools & compiler ###

RUN yum -y upgrade \
  && yum install -y centos-release-scl \
  scl-utils \
  which \
  bzip2 \
  make \
  git \
  wget \
  autoconf \
  autoconf-archive \
  automake \
  libtool \
  bash-completion \
  zlib-devel \
  && yum install -y $CASM_DEVTOOLSET-gcc \
  $CASM_DEVTOOLSET-gcc-c++ \
  $CASM_DEVTOOLSET-gcc-gfortran \
  && yum clean all -y
  
# Effects of `scl enable $CASM_DEVTOOLSET bash`
ENV MANPATH=/opt/rh/$CASM_DEVTOOLSET/root/usr/share/man:
ENV PERL5LIB=/opt/rh/$CASM_DEVTOOLSET/root/usr/lib64/perl5/vendor_perl:/opt/rh/$CASM_DEVTOOLSET/root/usr/lib/perl5:/opt/rh/$CASM_DEVTOOLSET/root/usr/share/perl5/vendor_perl
ENV X_SCLS=$CASM_DEVTOOLSET 
ENV PCP_DIR=/opt/rh/$CASM_DEVTOOLSET/root
ENV INFOPATH=/opt/rh/$CASM_DEVTOOLSET/root/usr/share/info
ENV PATH=/opt/rh/$CASM_DEVTOOLSET/root/usr/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/rh/$CASM_DEVTOOLSET/root/usr/lib64:/opt/rh/$CASM_DEVTOOLSET/root/usr/lib


### Stage 2: Install conda ###

# create casmgroup and casmuser
RUN groupadd casmgroup && useradd -g casmgroup casmuser

# create CASM_CONDA_DIR for casmuser
RUN mkdir -p /tmp \
  && mkdir -p $CASM_CONDA_DIR \
  && chown -R casmuser $CASM_CONDA_DIR

# install miniconda
USER casmuser
RUN curl -sSL https://repo.continuum.io/miniconda/Miniconda${CASM_PYTHON_VERSION:0:1}-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
  && mkdir -p $CASM_CONDA_DIR \
  && bash /tmp/miniconda.sh -bfp $CASM_CONDA_DIR \
  && PATH="$CASM_CONDA_DIR/bin:$PATH" \
  && rm -rf /tmp/miniconda.sh \
  && conda install -y \
    "python =$CASM_PYTHON_VERSION" \
    conda-build \
    anaconda-client \
  && conda update --all \
  && conda clean --all --yes

# include conda on path
ENV PATH=$CASM_CONDA_DIR/bin:$PATH

CMD ["/bin/bash"]
