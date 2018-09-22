ARG CASM_BRANCH=err
FROM bpuchala/casm-base:"$CASM_BRANCH"_condagcc

LABEL description="CASM, using Anaconda compiler tools"

# Make it easier to customize build by allowing build arguments
ARG CASM_CONDA_VERSION
ARG CASM_PYTHON_VERSION
ARG CASM_BOOST_VERSION
ARG CASM_BOOST_CONDAGCC_BUILD_STR
ARG CASM_CONDA_CHANNEL

ENV CASM_ENV_NAME=casm_"$CASM_CONDA_VERSION"_py"$CASM_PYTHON_VERSION"

USER casmuser

### Create conda environment for testing CASM builds ###

# create conda environment for building CASM
RUN conda create -n $CASM_ENV_NAME -y \
  -c $CASM_CONDA_CHANNEL -c defaults -c conda-forge -c prisms-center  \
  "python =$CASM_PYTHON_VERSION" \
  "casm-boost $CASM_BOOST_VERSION $CASM_BOOST_CONDAGCC_BUILD_STR"  \
  && conda clean -y -a

WORKDIR /home/casmuser

CMD ["/bin/bash"]
