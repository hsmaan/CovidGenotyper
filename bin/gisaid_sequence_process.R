library(ape)
library(Biostrings)
library(stringr)
library(DECIPHER)
library(data.table)
library(tidyverse)

# Load fasta file

setwd("../data")

all_files <- list.files()

gisaid_fastas <- grep("gisaid_cov2020_sequences", all_files, value = TRUE)

all_fastas <- readDNAStringSet(gisaid_fastas)

names(all_fastas) <- (str_split_fixed(names(all_fastas), fixed("|"), 3))[,2]

# Load reference fasta

ref_fasta <- readDNAStringSet("ncov_ref_NC_045512.fasta")

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
file.remove(grep("dec_aligned_plus_ref_filtered_", all_files, value = TRUE))
file.remove(grep("signal_aligned", all_files, value = TRUE))

# Remove sequences with high percentage of ambiguous nucleotides

ambg_freq_all <- which(((apply((letterFrequency(all_fastas, c("N", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V"))), 1, sum))/29000) > 0.001)

all_fastas <- all_fastas[-ambg_freq_all]

# Subset all files by available metadata and sample 10000 sequences 

all_fastas <- all_fastas[which(names(all_fastas) %in% meta_accession[,1])]
meta_accession <- meta_accession[which(meta_accession[,1] %in% names(all_fastas)),]
meta_accession <- data.frame("Accession" = as.character(meta_accession))
accession_order <- base::match(names(all_fastas), meta_accession[,1])
all_fastas <- all_fastas[order(accession_order)]
all_10000 <- sample(names(all_fastas), 10000)
all_fastas <- all_fastas[all_10000]

# Remove gaps and align

all_fastas <- RemoveGaps(all_fastas, removeGaps = "all", processors = NULL)
fasta_align <- AlignSeqs(all_fastas, iterations = 0, refinements = 0, processors = NULL)

# Subset and exclude UTRs and gaps 

fasta_mat <- as.matrix(fasta_align)
fasta_bin <- as.DNAbin(fasta_mat)
fasta_ungapped <- del.colgapsonly(fasta_bin, threshold = 0.95)
fasta_string <- fasta_ungapped %>% as.list %>% as.character %>% lapply(., paste0, collapse = "") %>% unlist %>% DNAStringSet
fasta_final <- subseq(fasta_string, start = 265, end = 29674)

writeXStringSet(fasta_string, file = paste("dec_aligned_fasta_filtered_", as.character(Sys.Date()), ".fasta", sep = ""))
writeXStringSet(fasta_string, file = "signal_aligned.fasta")

# Mask sites 

mask_sites <- c(187, 1059, 2094, 3037, 3130, 6990, 8022, 10323, 10741, 11074, 13408, 14786, 19684, 20148, 21137, 24034, 24378, 25563, 26144, 26461, 26681, 28077, 28826, 28854, 29700, 4050, 13402, 11083, 15324, 21575)
align_mat <- as.matrix(fasta_string)
align_mat_sub <- align_mat[, -mask_sites]
align_mat_bin <- as.DNAbin(align_mat_sub)
align_masked <- align_mat_bin %>% as.list %>% as.character %>% lapply(., paste0, collapse = "") %>% unlist %>% DNAStringSet
fasta_masked <- subseq(align_masked, start = 265, end = 29674)

# Get distance

fasta_dist <- dist.dna(as.DNAbin(fasta_masked), model = "K80", as.matrix = TRUE, pairwise.deletion = FALSE)

# Append ref to alignment, subset for seqs not equal to 29903 in length, and output 

fasta_align_ref <- c(ref_fasta, fasta_string)
fasta_align_ref <- fasta_align_ref[which(length(fasta_align_ref) == 29903)]
writeXStringSet(fasta_align_ref, file = paste("dec_aligned_plus_ref_filtered_", as.character(Sys.Date()), ".fasta", sep = ""))
pre_meta_sub <- pre_meta[which(pre_meta$Accession %in% names(fasta_align_ref)),]

# Output all objects 

save(fasta_align, file = paste("dec_unfiltered_aligned_fastas_", as.character(Sys.Date()), ".RData", sep = ""))
save(fasta_dist, file = paste("dec_fasta_dist_", as.character(Sys.Date()), ".RData", sep = ""))
save(fasta_string, file = paste("dec_aligned_filtered_", as.character(Sys.Date()), ".RData", sep = ""))
save(pre_meta_sub, file = paste("covid_filtered_meta_", as.character(Sys.Date()), ".RData", sep = ""))

