FROM rocker/tidyverse:4.2
MAINTAINER rokita@chop.edu
WORKDIR /rocker-build/

RUN RSPM="https://packagemanager.rstudio.com/cran/2022-10-07" \
  && echo "options(repos = c(CRAN='$RSPM'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site
  
COPY scripts/install_bioc.r .

COPY scripts/install_github.r .

### Install apt-getable packages to start
#########################################
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog

# Add curl, bzip2 and some dev libs
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    curl \
    bzip2 \
    zlib1g \
    libbz2-dev \
    liblzma-dev \
    libreadline-dev

# libmagick++-dev is needed for coloblindr to install
RUN apt-get -y --no-install-recommends install \
    libgdal-dev \
    libudunits2-dev \
    libmagick++-dev

# Required for installing pdftools, which is a dependency of gridGraphics
RUN apt-get -y --no-install-recommends install \
    libpoppler-cpp-dev

# Install java
RUN apt-get update && apt-get -y --no-install-recommends install \
   default-jdk \
   libxt6
   
# cmakeis needed for ggpubr to install
RUN apt-get -y --no-install-recommends install \
    cmake

# install R packages from CRAN
RUN install2.r \
	BiocManager \
  ComplexHeatmap \
  data.table \
  ggpubr \
  ggthemes \
  msigdbr \
	optparse \
	pheatmap \
	RColorBrewer \
	survival \
  survMisc \
  survminer \
  tidytext

	
# install R packages from GitHub

# Maftools
RUN ./install_github.r \
	PoisonAlien/maftools

# Patchwork for plot compositions
RUN ./install_github.r  'thomasp85/patchwork' --ref 'c67c6603ba59dd46899f17197f9858bc5672e9f4'



# Install pip3 and python reqs for oncokb
RUN apt-get -y --no-install-recommends install \
    python3-pip python3-dev
RUN pip3 install \
  "matplotlib==3.1.2" \
  "kiwisolver==1.2.0" \
  "requests==2.27.1" \
  "urllib3==1.26.8"

# Install oncokb
RUN git clone https://github.com/oncokb/oncokb-annotator.git /home/oncokb-annotator


WORKDIR /rocker-build/

ADD Dockerfile .
