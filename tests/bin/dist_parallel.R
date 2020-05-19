library(parallel)
library(Biostrings)
library(ggplot2)
library(ggthemes)
library(ape)

# All cores

cores <- detectCores()

# Load aligned sequences

setwd("../../data")

file_list <- list.files()
gisaid_fastas <- grep("dec_aligned_fasta_filtered", file_list, value = TRUE)
aligned <- readDNAStringSet(gisaid_fastas)

setwd("../tests/bin")

# 

aligned_1 <- aligned[1]
aligned_M1 <- aligned[-1]

dna_dist <- function(d1, d2) {
  
  d1 <- DNAStringSet(d1)
  d2 <- DNAStringSet(d2)
  d3 <- c(d1, d2)
  d3_dist <- dist.dna(d3)
  return(d3_dist)
  
}

test <- lapply(aligned_M1, function(x) dna_dist(aligned_1, x))

lapply(aligned_M1, print)
