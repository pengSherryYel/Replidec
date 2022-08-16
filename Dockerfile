FROM mambaorg/micromamba:0.22.0

LABEL author="Jinlong Ru"

USER root

# Install packages into `base` env, because nextflow can't use other envs.
RUN micromamba install -n base -c denglab -c conda-forge -c bioconda replidec && \
    micromamba clean --all --yes

# Activate conda env during docker build
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH

RUN cd /opt/conda/lib/python3.10/site-packages/Replidec && \
    wget https://zenodo.org/record/6975142/files/db_v0.2.3.tar.gz && \
    tar xzf db_v0.2.3.tar.gz && \
    rm db_v0.2.3.tar.gz

