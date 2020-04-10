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
library(tidyverse)

# Load fasta file

setwd("../data")

all_files <- list.files()

gisaid_fastas <- grep("gisaid_cov2020_sequences", all_files, value = TRUE)

all_fastas <- readDNAStringSet(gisaid_fastas)

names(all_fastas) <- (str_split_fixed(names(all_fastas), fixed("|"), 3))[,2]

# Load metadata 

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

pre_meta <- loadRData(grep("covid_meta", all_files, value = TRUE))

meta_accession <- data.frame("Accession" = as.character(pre_meta$Accession))

# Remove old files

file.remove(grep("dec_aligned_fasta_filtered_", all_files, value = TRUE))
file.remove(grep("dec_unfiltered_aligned_fastas_", all_files, value = TRUE))
file.remove(grep("dec_fasta_dist_", all_files, value = TRUE))
file.remove(grep("dec_aligned_filtered_", all_files, value = TRUE))
file.remove(grep("covid_filtered_meta_", all_files, value = TRUE))

# Subset all files by available metadata

all_fastas <- all_fastas[which(names(all_fastas) %in% meta_accession[,1])]

meta_accession <- meta_accession[which(meta_accession[,1] %in% names(all_fastas)),]

meta_accession <- data.frame("Accession" = as.character(meta_accession))

accession_order <- base::match(names(all_fastas), meta_accession[,1])

all_fastas <- all_fastas[order(accession_order)]

pre_meta_sub <- pre_meta[which(pre_meta$Accession %in% names(all_fastas)),]

# Remove gaps

all_fastas <- RemoveGaps(all_fastas, removeGaps = "all", processors = NULL)

# Align 

fasta_align <- AlignSeqs(all_fastas, iterations = 0, refinements = 0)

# Subset and exclude UTRs and gaps 

fasta_mat <- as.matrix(fasta_align)
fasta_bin <- as.DNAbin(fasta_mat)
fasta_ungapped <- del.colgapsonly(fasta_bin, threshold = 0.95)
fasta_string <- fasta_ungapped %>% as.list %>% as.character %>% lapply(., paste0, collapse = "") %>% unlist %>% DNAStringSet

fasta_final <- subseq(fasta_string, start = 265, end = 29674)

writeXStringSet(fasta_string, file = paste("dec_aligned_fasta_filtered_", as.character(Sys.Date()), ".fasta", sep = ""))

# Get distance

fasta_dist <- dist.dna(as.DNAbin(fasta_final), model = "K80", as.matrix = TRUE, pairwise.deletion = FALSE)

# Output both objects 

save(fasta_align, file = paste("dec_unfiltered_aligned_fastas_", as.character(Sys.Date()), ".RData", sep = ""))
save(fasta_dist, file = paste("dec_fasta_dist_", as.character(Sys.Date()), ".RData", sep = ""))
save(fasta_string, file = paste("dec_aligned_filtered_", as.character(Sys.Date()), ".RData", sep = ""))
save(pre_meta_sub, file = paste("covid_filtered_meta_", as.character(Sys.Date()), ".RData", sep = ""))

