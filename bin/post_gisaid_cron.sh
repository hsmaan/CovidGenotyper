#!/bin/bash

cd /home/hmaan/Documents/covid_shiny/bin

# Run processing scripts sequentially

sh snp_sites_process.sh

Rscript maf_sites_out.R

Rscript preprocess_plot_data.R

# Delete old docker image and rebuild new. Push new image to google container reg

docker image rm gcr.io/cgt-100/cgt 

docker system prune --force

cd ..

docker build -t gcr.io/cgt-100/cgt -f app.Dockerfile .

docker push gcr.io/cgt-100/cgt:latest

# Reinstantiate vm instance

cd bin

Rscript gce_app_deploy.R 

echo "Update completed"





