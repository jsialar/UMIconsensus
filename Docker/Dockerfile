FROM ubuntu:24.04

RUN apt-get update && \
    apt-get install -y wget && \
    apt-get install -y gcc && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py312_25.7.0-2-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda
ENV PATH=/opt/conda/bin:${PATH}

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda update -y conda

COPY environment.yml .

RUN conda env update -n root -f environment.yml && \
    conda clean --all

WORKDIR "/opt"
