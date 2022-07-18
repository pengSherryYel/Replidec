FROM mambaorg/micromamba:0.22.0

LABEL author="Jinlong Ru"

USER root

# Install packages into `base` env, because nextflow can't use other envs.
RUN micromamba install -n base -c denglab -c conda-forge -c bioconda replidec && \
    micromamba clean --all --yes

# Activate conda env during docker build
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH
