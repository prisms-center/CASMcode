ARG CASM_CONDA_VERSION:err
FROM bpuchala/casm-build-condagcc:$CASM_CONDA_VERSION

LABEL description="CASM, using Anaconda compiler tools"

# Make it easier to customize build by allowing build arguments
ARG CASM_CONDA_VERSION

ENV CASM_ENV_NAME=casm_"$CASM_CONDA_VERSION"_py"$CASM_PYTON_VERSION"

USER casmuser

### Create conda environment for CASM ###

# create conda environment for CASM
RUN conda create -n $CASM_ENV_NAME -y \
  python=$CASM_PYTHON_VERSION \
  casm-condagcc=$CASM_CONDA_VERSION \
  && conda clean -y -a

WORKDIR /home/casmuser

CMD ["/bin/bash"]
