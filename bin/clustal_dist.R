library(ape)
library(adephylo)
library(ade4)
library(phylobase)
library(geiger)
library(phytools)
library(muscle)
library(Biostrings)
library(stringr)
library(DECIPHER)
library(phytools)
library(stringi)
library(data.table)

dist_get <- function(align) {
  
  # dec_align_pre <- dec_align[1:length(covid_align)]
  # dec_align_post <- dec_align[length(covid_align):length(dec_align)]
  # dist_get <- (seq1, seq2) {
  #   
  #   dist <- dist.dna(as.DNAbin())
  # }
  dec_dist <- dist.dna(as.DNAbin(align), model = "K80", as.matrix = TRUE, pairwise.deletion = FALSE)
  colnames(dec_dist) <- (str_split_fixed(colnames(dec_dist), fixed("."), 2)[,1])
  rownames(dec_dist) <- (str_split_fixed(rownames(dec_dist), fixed("."), 2)[,1])
  return(dec_dist)
    
}

