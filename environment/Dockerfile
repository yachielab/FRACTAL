# hash:sha256:0c19494bcb0f18971ad3a94779994f89f6c6c227b833f49e69e1800398eca5fb
FROM registry.codeocean.com/codeocean/miniconda3:4.7.10-python3.7-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN conda install -y \
        biopython==1.76 \
    && conda clean -ya

COPY postInstall /
RUN /postInstall

ENV PATH /root/capsule/code/FRACTAL:/root/capsule/code/FRACTAL/example/script:$PATH