FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive
ARG USER_ID=root
ENV USER_ID $USER_ID

RUN  apt-get update \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/*

USER $USER_ID
WORKDIR /home/$USER_ID
ENV HOME /home/$USER_ID

## Git install
RUN apt-get -y update
RUN apt-get -y install git

RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base
## enable R package install
RUN chmod a+w /usr/local/lib/R/site-library

RUN apt-get update && apt-get upgrade -y \
    && apt-get install -y build-essential cmake git \
       ssh g++ r-cran-dplyr r-cran-crayon r-cran-rcpp \
       r-cran-devtools r-cran-ggplot2 r-cran-ggraph \
       r-cran-knitr r-cran-r.utils r-cran-tidygraph \
       r-cran-hexbin r-cran-rcolorbrewer r-cran-cli \
    && apt-get clean

RUN R -e 'install.packages(c("ggmuller"))'
RUN R -e 'devtools::install_github("caravagnalab/ProCESS", ref="label_tour")'
RUN R -e 'install.packages("tidyverse")'
RUN R -e 'devtools::install_github("thomasp85/patchwork")'

RUN apt-get install -y libncurses-dev
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y make

## install samtools
RUN wget 'https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2'
RUN tar -xvjf samtools-1.20.tar.bz2
WORKDIR samtools-1.20/
RUN sh ./configure --prefix=/home --disable-bz2 --disable-lzma
RUN make
RUN make install
ENV PATH="/home/bin:$PATH"

# install python
##RUN apt-get update && apt-get install -y --no-install-recommends build-essential python3.6 python3-pip python3-setuptools python3-dev
RUN apt-get -y install python3 python3-pip python3-setuptools python3-dev
