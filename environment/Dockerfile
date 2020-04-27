# hash:sha256:31bfe4d8e1e914f67bd995af2709190aa47f03c5b9e33bd9f0e740b34524d5f5
FROM registry.codeocean.com/codeocean/miniconda3:4.7.10-python3.7-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN conda install -y --channel=bioconda \
        biopython==1.76 \
        epa-ng==0.3.5 \
        rapidnj==v2.3.2 \
    && conda clean -ya

COPY postInstall /
RUN /postInstall
