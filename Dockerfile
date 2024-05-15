FROM rocker/tidyverse:4.4.0
MAINTAINER rokita@chop.edu
WORKDIR /rocker-build/

#RUN RSPM="https://packagemanager.rstudio.com/cran/2022-10-07" \
#  && echo "options(repos = c(CRAN='$RSPM'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site
  
#COPY scripts/install_bioc.r .

#COPY scripts/install_github.r .

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

# Set the Bioconductor repository as the primary repository
RUN R -e "options(repos = BiocManager::repositories())"

# Install BiocManager and the desired version of Bioconductor
RUN R -e "install.packages('BiocManager', dependencies=TRUE)"
RUN R -e "BiocManager::install(version = '3.19')"

# Install packages
RUN R -e 'BiocManager::install(c( \
	"BiocManager", \
  "data.table", \
  "ggpubr", \
  "ggthemes", \
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
  "GSVA" \
))'

# install R packages from CRAN
#RUN install2.r \
#	BiocManager \
#  data.table \
#  ggpubr \
#  ggthemes \
#  msigdbr \
#  openxlsx \
#	optparse \
#	pheatmap \
#	RColorBrewer \
#	survival \
# survMisc \
#  survminer \
#  tidytext
  
  
# install R packages from Bioconductor 
#RUN ./install_bioc.r \
#  ComplexHeatmap \
#  GSVA 

# install R packages from GitHub

# Maftools
#RUN ./install_github.r \
#	PoisonAlien/maftools


# Patchwork for plot compositions
RUN R -e "remotes::install_github('thomasp85/patchwork', ref = 'c67c6603ba59dd46899f17197f9858bc5672e9f4', dependencies = TRUE)"
RUN R -e "remotes::install_github('clauswilke/colorblindr', ref = '90d64f8fc50bee7060be577f180ae019a9bbbb84', dependencies = TRUE)"
RUN R -e "remotes::install_github('coolbutuseless/ggpattern', ref = 'bc9e4cd1271674a537bf3406663765648e3963bd', dependencies = TRUE)"
RUN R -e "remotes::install_github('PoisonAlien/maftools', ref = 'ecaf525b95449b719fa62c1e693aa67a9356b344', dependencies = TRUE)"

# Install pip3 and python reqs for oncokb
RUN apt-get update
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
