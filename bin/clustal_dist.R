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

dist_get <- function(fasta, align) {
  
  covid_seq <- readDNAStringSet(fasta)
  covid_align <- align
  dec_align <- AlignProfiles(covid_seq, covid_align, restrict = c(-100, 2, 10))
  dec_dist <- dist.dna(as.DNAbin(dec_align), model = "K80", as.matrix = TRUE, pairwise.deletion = FALSE)
  return(dec_dist)
    
}

