# hash:sha256:aaa44c3ffa6ddf75b78ccc025219f559b7bd9c1c5338b5048421fccea709fee2
FROM registry.codeocean.com/codeocean/miniconda3:4.8.2-python3.8-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN conda install -y --channel=bioconda \
        python==3.7.1 \
        biopython==1.76 \
        epa-ng==0.3.5 \
        hmmer==3.3 \
        mafft==7.464 \
        rapidnj==v2.3.2 \
        seqkit==0.11.0 \
        trimal==1.4.1 \
    && conda clean -ya

COPY postInstall /
RUN /postInstall
