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

dist_get <- function(fasta, align, seqnames) {
  
  covid_seq <- readDNAStringSet(fasta)
  covid_seq <- as.DNAbin(covid_seq)
  names(align) <- seqnames
  covid_align <- as.DNAbin(align)
  clustal_align <- AlignProfiles(x = covid_seq, y = covid_align, pw.gapopen = 10, pw.gapext = 0.1, gapopen = 10, gapext = 0.2, exec = "/home/hmaan/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2")
  clustal_dist <- dist.dna(clustal_align, model = "K80", as.matrix = TRUE, pairwise.deletion = FALSE)
  return(clustal_dist)
    
}

dist_get_null <- function(align) {
  
  covid_align <- as.DNAbin(align)
  covid_dist <- dist.dna(covid_align, model = "K80", as.matrix = TRUE, pairwise.deletion = FALSE)
  return(covid_dist)
  
}

test <- readDNAStringSet("data/covid_ncbi_repeat.fasta")
test_align <- AlignSeqs(test, iterations = 0, refinements = 0)

test2 <- readDNAStringSet("data/ncbi_test_seq.fasta")
  
clustal_dist <- dist.dna(as.DNAbin(test_align), model = "K80", as.matrix = TRUE, pairwise.deletion = FALSE)
u <- umap_process(as.numeric.matrix(clustal_dist))
plot(u)

test_3 <- AlignProfiles(test2, test_align, restrict = c(-100, 2, 10))
dist_plu
