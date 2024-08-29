FROM rocker/tidyverse:4.4.0
MAINTAINER rokita@chop.edu
WORKDIR /rocker-build/

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
   
# cmake is needed for ggpubr to install
RUN apt-get -y --no-install-recommends install \
    cmake

# add bedtools
RUN apt-get install -y bedtools && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the Bioconductor repository as the primary repository
RUN R -e "options(repos = BiocManager::repositories())"

# Install BiocManager and the desired version of Bioconductor
RUN R -e "install.packages('BiocManager', dependencies=TRUE)"
RUN R -e "BiocManager::install(version = '3.19')"

# Install packages
RUN R -e 'BiocManager::install(c( \
	"BiocManager", \
	"biomaRt", \
  "data.table", \
  "ggpubr", \
  "ggthemes", \
  "maftools", \
  "msigdbr", \
  "openxlsx", \
	"optparse", \
	"pheatmap", \
	"RColorBrewer", \
	"survival", \
  "survMisc", \
  "survminer", \
  "tidytext", \
  "ComplexHeatmap", \
  "GSVA", \
  "R.utils" \
))'

# install R packages from GitHub

RUN R -e "remotes::install_github('thomasp85/patchwork', ref = '1cb732b129ed6a65774796dc1f618558c7498b66', dependencies = TRUE)"
RUN R -e "remotes::install_github('clauswilke/colorblindr', ref = '90d64f8fc50bee7060be577f180ae019a9bbbb84', dependencies = TRUE)"
RUN R -e "remotes::install_github('coolbutuseless/ggpattern', ref = 'bc9e4cd1271674a537bf3406663765648e3963bd', dependencies = TRUE)"

# Install pip3 and python reqs for oncokb
RUN apt-get update
RUN apt-get -y --no-install-recommends install \
    python3-pip python3-dev
RUN python3 -m pip install --upgrade pip
RUN pip3 install \
  "matplotlib==3.1.2" \
  "kiwisolver==1.2.0" \
  "requests==2.27.1" \
  "urllib3==1.26.8"

# Install oncokb
RUN git clone https://github.com/oncokb/oncokb-annotator.git /home/oncokb-annotator


WORKDIR /rocker-build/

ADD Dockerfile .
