#!/bin/bash

cd /home/hmaan/Documents/covid_shiny/bin

# Run processing scripts sequentially

sh snp_sites_process.sh

Rscript maf_sites_out.R

Rscript preprocess_plot_data.R

