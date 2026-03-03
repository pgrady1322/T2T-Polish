# syntax=docker/dockerfile:1
FROM --platform=linux/amd64 ubuntu:22.04

# Avoid tzdata interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# 1. System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential ca-certificates wget git cmake \
    python3 python3-pip python3-dev \
    curl unzip default-jre \
    zlib1g-dev libbz2-dev liblzma-dev \
    libcurl4-gnutls-dev libssl-dev \
    samtools bcftools && \
    rm -rf /var/lib/apt/lists/*

# 2. Install Miniconda
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/conda.sh && \
    bash /tmp/conda.sh -b -p $CONDA_DIR && rm /tmp/conda.sh
ENV PATH=$CONDA_DIR/bin:$PATH

# 3. Create and populate the bioinformatics environment
RUN conda update -n base -c defaults conda && \
    conda create -y -n apv4-env -c defaults -c conda-forge -c bioconda \
      python=3.12 \
      meryl \
      winnowmap \
      kmer-jellyfish \
      genomescope2 \
      samtools \
      pysam \
      tqdm \
      pb-falconc && \
    conda clean --all -y

# 4. Build MerFin v1.1 from source
RUN git clone --depth 1 --branch v1.1 --recurse-submodules \
        https://github.com/arangrhie/merfin.git /opt/Merfin && \
    git -C /opt/Merfin submodule update --init --recursive && \
    cd /opt/Merfin/src && make -j$(nproc) && \
    install -m 0755 /opt/Merfin/build/bin/merfin $CONDA_DIR/envs/apv4-env/bin/merfin && \
    rm -rf /opt/Merfin

# 5. Install Merqury from source
RUN git clone --depth 1 https://github.com/marbl/merqury.git /opt/Merqury && \
    install -m 0755 /opt/Merqury/merqury.sh $CONDA_DIR/envs/apv4-env/bin/merqury.sh && \
    ln -sf $CONDA_DIR/envs/apv4-env/bin/merqury.sh $CONDA_DIR/envs/apv4-env/bin/merqury && \
    rm -rf /opt/Merqury

# 6. Copy the T2T-Polish package and install it
WORKDIR /opt/t2t-polish
COPY pyproject.toml .
COPY t2t_polish/ t2t_polish/

# Install in the conda env so the 't2t-polish' console script is available
RUN $CONDA_DIR/envs/apv4-env/bin/pip install --no-deps .

# 7. Ensure the conda env is activated by default
SHELL ["/bin/bash", "-c"]
RUN echo "conda activate apv4-env" >> ~/.bashrc
ENV PATH=$CONDA_DIR/envs/apv4-env/bin:$PATH

# 8. Entrypoint
ENTRYPOINT ["t2t-polish"]