FROM centos:7

LABEL description="CASM build & run container, using Anaconda compiler tools"

# Make it easier to customize build by allowing build arguments
ARG PYTHON_VERSION
ARG CONDAGCC_VERSION

ENV CASM_PYTHON_VERSION=$PYTHON_VERSION
ENV CASM_CONDAGCC_VERSION=$CONDAGCC_VERSION
ENV CASM_CONDA_DIR=/home/casmuser/.local/conda

### yum installs ###

RUN yum -y upgrade \
  && yum install -y which \
  bzip2 \
  bash-completion \
  make \
  git \
  wget \
  autoconf \
  automake \
  libtool \
  && yum clean all -y

### Install conda ###

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
    "gcc_linux-64 =$CASM_CONDAGCC_VERSION" \
    "gxx_linux-64 =$CASM_CONDAGCC_VERSION" \
    "gfortran_linux-64 =$CASM_CONDAGCC_VERSION"\
  && conda update --all \
  && conda clean --all --yes

# include conda on path
ENV PATH=$CASM_CONDA_DIR/bin:$PATH

WORKDIR /home/casmuser

CMD ["/bin/bash"]
