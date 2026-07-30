# Start from the official, rock-solid R 4.5.2 base image
FROM rocker/r-ver:4.5.2

# 1. Install critical system-level libraries required by genomics R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    libglpk-dev \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libgit2-dev \
    libgsl-dev \
    libperl-dev \
    make \
    cmake \
    gfortran \
    automake \
    autoconf \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# 2. Install BiocManager first to manage Bioconductor components
RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

# 3. Install core CRAN and Bioconductor packages with all their dependencies
RUN Rscript -e "BiocManager::install(c( \
    'tidyverse', \
    'ggplot2', \
    'GenomicRanges' \
), update=FALSE, ask=FALSE)"

RUN Rscript -e "BiocManager::install(c( \
    'R.utils', \
    'data.table', \
    'optparse', \
    'GenomeInfoDb', \
    'IRanges', \
    'BSgenome.Hsapiens.UCSC.hg38' \
), update=FALSE, ask=FALSE)"

RUN Rscript -e "install.packages('parallel')"

#  Clean up temporary directory file spaces to keep container layers lightweight
WORKDIR /
RUN rm -rf /tmp/build
