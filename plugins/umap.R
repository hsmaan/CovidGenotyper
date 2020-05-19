library(Biostrings)
library(stringr)
library(data.table)
library(uwot)

# Read args

args <- commandArgs(trailingOnly = TRUE)

fasta_file <- as.character(args[1])
align_file <- as.character(args[2])
dist_file <- as.character(args[3])
meta_file <- as.character(args[4])

# Load files and format 

con_fasta <- readDNAStringSet(fasta_file)

align_fasta <- readDNAStringSet(align_file)

dist_dt <- fread(dist_file, stringsAsFactors = FALSE)
dist_mat <- as.matrix(dist_dt[, -1, with = FALSE])
rm(dist_dt)
gc()
colnames(dist_mat) <- colnames(dist_dt[, -1])
rownames(dist_mat) <- as.vector(dist_dt[[1]])

meta_dt <- fread(meta_file, stringsAsFactors = FALSE)
meta_df <- as.data.frame(meta_dt)
meta_df_sub <- meta_df[,c("gisaid_epi_isl", "region", "country", "country_exposure", "date")]
colnames(meta_df_sub) <- c("Accession", "Region", "Country", "Travel_history", "Date")
meta_df_sub$Travel_history <- ifelse(meta_df_sub$Country == meta_df_sub$Travel_history, "Not available", meta_df_sub$Travel_history)

# 

split_dist <- function(d1, d2) {
  
  d1 <- DNAStringSet(d1)
  d2 <- DNAStringSet(d2)
  dna_bound <- c(d1, d2)
  bound_dist <- dist_get(dna_bound)
  bound_df <- as.data.frame(bound_dist[2])
  colnames(bound_df) <- names(dna_bound)[1]
  rownames(bound_df) <- names(dna_bound)[2]
  return(bound_df)
  
}

align_bind <- function(full_dist, row_list) {
  
  new_row <- row_list[[1]]
  name <- row_list[[2]]
  new_col <- c(new_row, 0)
  dist_new <- rbind(full_dist, new_row)
  rownames(dist_new)[nrow(dist_new)] <- name
  dist_new <- cbind(dist_new, new_col)
  colnames(dist_new)[ncol(dist_new)] <- name
  return(dist_new)
}
