FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for IKMB HLA imputation pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/hla-pipe-1.0/bin:$PATH
