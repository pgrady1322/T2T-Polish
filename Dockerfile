# syntax=docker/dockerfile:1
FROM --platform=linux/amd64 ubuntu:20.04

# Avoid tzdata interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# 1. System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential ca-certificates wget git cmake python3 python3-pip \
    curl unzip default-jre zlib1g-dev libbz2-dev liblzma-dev \
    samtools bcftools && \
    rm -rf /var/lib/apt/lists/*

 RUN apt-get update && apt-get install -y \
    python3-pip python3-dev \
    libbz2-dev liblzma-dev zlib1g-dev \
    libcurl4-gnutls-dev libssl-dev && \
    rm -rf /var/lib/apt/lists/*

# Install Python requirements for APv4.py
RUN pip3 install pysam tqdm

# 2. Install Miniconda
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/conda.sh && \
    bash /tmp/conda.sh -b -p $CONDA_DIR && rm /tmp/conda.sh
ENV PATH=$CONDA_DIR/bin:$PATH

# 3. Create and populate the bioinformatics environment
RUN conda update -n base -c defaults conda && \
    conda create -y -n apv4-env -c defaults -c conda-forge -c bioconda \
      python=3.6 \
      deepvariant=1.4.0 \
      meryl \
      winnowmap \
      jellyfish \
      genomescope2 \
      samtools \
      pb-falconc && \
    conda clean --all -y

# 4. Build MerFin v1.1 from source
RUN git clone https://github.com/arangrhie/merfin.git /opt/Merfin && \
    cd /opt/Merfin && git checkout v1.1 && \
    cd src && make -j$(nproc) && \
    cp /opt/Merfin/build/bin/merfin $CONDA_DIR/envs/apv4-env/bin/

# 5. Install Merqury from source
RUN git clone https://github.com/marbl/merqury.git /opt/Merqury && \
    cd /opt/Merqury && \
    mkdir -p /opt/Merqury/bin && \
    cp merqury.sh /opt/Merqury/bin/merqury && \
    chmod +x /opt/Merqury/bin/merqury && \
    ln -s /opt/Merqury/bin/merqury $CONDA_DIR/envs/apv4-env/bin/merqury

# 7. Copy your APv4 script into PATH
COPY APv4.py /usr/local/bin/APv4.py
RUN chmod +x /usr/local/bin/APv4.py

# 8. Ensure the conda env is activated by default
SHELL ["/bin/bash", "-lc"]
RUN echo "conda activate apv4-env" >> ~/.bashrc

# 9. Optional default entrypoint
ENTRYPOINT ["bash", "-lc", "APv4.py"]