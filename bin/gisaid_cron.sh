#!/bin/bash

cd /home/hmaan/covid_shiny/bin

# Run metadata and sequence processing scripts

Rscript --verbose metadata_process.R

Rscript --verbose gisaid_sequence_process.R

# Run snp processing scripts

sh snp_sites_process.sh

Rscript --verbose maf_sites_out.R

# Preprocess the plot data

Rscript --verbose preprocess_plot_data.R

# Delete old docker image and rebuild new. Push new image to google container reg

docker image rm gcr.io/bo-wang-274419/cgt 

docker system prune --force

cd ..

docker build -t gcr.io/bo-wang-274419/cgt -f app.Dockerfile .

docker push gcr.io/bo-wang-274419/cgt

