#!/bin/bash

module load R/3.6.1

sbatch --mem=16G -J cv_meta -t 5-00:00:00 -c 12 Rscript --verbose metadata_process.R

sleep 1h 

sbatch -p himem --mem=50G -J cv_msa -t 5-00:00:00 -c 12 Rscript --verbose gisaid_sequence_process.R
