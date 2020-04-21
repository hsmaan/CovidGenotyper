library(data.table)
library(plyr)
library(RColorBrewer)
library(profvis)
library(parallel)

# Load testing data 

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

setwd("../data")

file_list <- list.files() 

pre_aligned_filtered <- loadRData(grep("dec_aligned_filtered", file_list, value = TRUE))
pre_meta <- loadRData(grep("covid_filtered_meta", file_list, value = TRUE))
pre_dist <- loadRData(grep("dec_fasta_dist", file_list, value = TRUE))
meta <- as.data.frame(pre_meta[,c("Accession", "Region", "Geo_Location", "Datetime")])
vars_freq <- loadRData(grep("var_freq_sub*", file_list, value = TRUE))

setwd("..")

# Source global

source("R/global.R")

setwd("data")

# DNA dist using ape

pre_aligned_dnabin <- as.DNAbin(pre_aligned_filtered)

test1_dist <- dist.dna(pre_aligned_dnabin, model = "K80", as.matrix = TRUE, pairwise.deletion = FALSE)