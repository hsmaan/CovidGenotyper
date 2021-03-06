# Rocker - R + Shiny

FROM rocker/shiny:3.6.3

# Set tmpdir

ENV TMPDIR=/home/hmaan/tmp

# General usage system libraries

RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev 

# Make directories 

RUN mkdir -p /bin/
RUN mkdir -p /data/
RUN mkdir -p /srv/shiny-server/R
RUN mkdir -p /srv/shiny-server/data

# Copy package install files

COPY ./bin/packages_install.R /bin/

# Install R packages

RUN Rscript /bin/packages_install.R


