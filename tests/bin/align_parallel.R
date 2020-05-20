library(parallel)
library(Biostrings)
library(DECIPHER)
library(ggplot2)
library(ggthemes)

# All cores

cores <- detectCores()

# Load raw alignment data

setwd("../../data")

file_list <- list.files()
gisaid_fastas <- grep("gisaid_cov2020_sequences", file_list, value = TRUE)
unaligned <- readDNAStringSet(gisaid_fastas)

setwd("../tests/bin")

# Define profile alignment function

align_profs <- function(x, y) {
  
  AlignProfiles(x, y, processors = NULL)
  
}

# Define full alignment function

align_full <- function(fastas_sub) {
  
  all_fastas <- RemoveGaps(fastas_sub, removeGaps = "all", processors = NULL)
  fasta_align <- AlignSeqs(all_fastas, iterations = 0, refinements = 0, processors = NULL)
  fasta_mat <- as.matrix(fasta_align)
  fasta_bin <- as.DNAbin(fasta_mat)
  fasta_ungapped <- del.colgapsonly(fasta_bin, threshold = 0.95)
  fasta_string <- fasta_ungapped %>% as.list %>% as.character %>% lapply(., paste0, collapse = "") %>% unlist %>% DNAStringSet
  return(fasta_string)
}
# Define full alignment test function

align_subsets <- function(fastas, div) {
  
  t1 <- Sys.time()
  ua_subsets <- split(fastas, ceiling(seq_along(aligned)/div))
  subsets_aligned <- mclapply(ua_subsets, align_full, mc.cores = cores)
  align_comp <- base::Reduce(align_profs, subsets_aligned)
  t2 <- Sys.time()
  tdiff <- t2 - t1
  return(list(align_comp, tdiff))
  
}

# Define chunk sizes to test

chunks <- c(10, 25, 50, 100, 250, 500)

# Lapply across all tests

align_results <- lapply(chunks, function(x) align_subsets(unaligned, x))

# Get results dataframe

align_length <- lapply(align_results[[1]], function(x) mean(width(x)))
times <- align_results[[2]]

align_df <- data.frame("chunk_size" = chunks, "time" = times, "align_size" = align_length)

# Plot results

setwd("../data")

ggplot(data = align_df, aes(x = chunk_size, y = time)) +
  theme_min() +
  geom_point() +
  labs(x = "Chunk size", y = "Time (minutes)")

ggsave("figures/align_test_size_time.pdf", device = "pdf", height = 7, width = 7, units = "in")

ggplot(data = align_df, aes(x = chunk_size, y = align_length)) +
  theme_min() +
  geom_point() +
  labs(x = "Chunk size", y = "Average sequence length")

ggsave("figures/align_test_size_time.pdf", device = "pdf", height = 7, width = 7, units = "in")


  



  
