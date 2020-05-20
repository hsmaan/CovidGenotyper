library(Biostrings)
library(stringr)
library(DECIPHER)
library(ape)
library(data.table)
library(ggplot2)
library(ggthemes)

# Functions from CovidGenotyper/R/global.R  

align_get <- function(fasta, align) {
  
  covid_seq <- fasta
  if (length(covid_seq) > 10) {
    stop("Too many sequences in fasta file, only up to 10 viral genomes allowed at a time.")
  }
  if (length(covid_seq) < 1) {
    stop("No sequences found in fasta file.")
  }
  if (length(which((mclapply(covid_seq, length, mc.cores = cores)) < 29000)) > 0) {
    stop("One or more sequences not complete (length < 29000 nucleotides).")
  }
  covid_align <- align
  covid_seq <- RemoveGaps(covid_seq, removeGaps = "all", processors = NULL)
  fasta_final <- AlignProfiles(covid_align, covid_seq, processors = NULL)
  return(fasta_final)
  
}

dist_get <- function(align) {
  
  mask_sites <- c(187, 1059, 2094, 3037, 3130, 6990, 8022, 10323, 10741, 11074, 13408, 14786, 19684, 20148, 21137, 24034, 24378, 25563, 26144, 26461, 26681, 28077, 28826, 28854, 29700, 4050, 13402, 11083, 15324, 21575)
  align_mat <- as.matrix(align)
  align_mat_sub <- align_mat[, -mask_sites]
  align_mat_bin <- as.DNAbin(align_mat_sub)
  align_masked <- align_mat_bin %>% as.list %>% as.character %>% lapply(., paste0, collapse = "") %>% unlist %>% DNAStringSet
  align_trim <- subseq(align_masked, start = 265, end = 29674)
  dec_dist <- dist.dna(as.DNAbin(align_trim), model = "K80", as.matrix = TRUE, pairwise.deletion = FALSE)
  colnames(dec_dist) <- (str_split_fixed(colnames(dec_dist), fixed("."), 2)[,1])
  rownames(dec_dist) <- (str_split_fixed(rownames(dec_dist), fixed("."), 2)[,1])
  return(dec_dist)
  
}

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

align_update <- function(dist_mat, min_id, name) {
  
  dist_min_row <- dist_mat[min_id[1],]
  return(list(dist_min_row, name))
  
}

dist_get_heur <- function(align, fasta, dist) {
  
  dist_ret0 <- dist
  new_length <- length(fasta)
  align_new <- align[(length(align) - new_length + 1):length(align)]
  align_minus <- align[1:(length(align) - new_length)]
  align_new_list <- as.list(align_new)
  
  split_test <- mclapply(align_new, function(x) mclapply(align_minus, function(y) split_dist(x, y), mc.cores = cores), mc.cores = cores)
  split_test_reduced <- mclapply(split_test, function(x) Reduce(cbind, x), mc.cores = cores)
  split_test_reduced_mins <- mclapply(split_test_reduced, min, mc.cores = cores)
  split_test_reduced_min_id <- mclapply(split_test_reduced, which.min, mc.cores = cores)
  split_test_reduced_names <- as.list(names(split_test_reduced_mins))
  
  if (max(unlist(split_test_reduced_mins)) > 1e-4) {
    dist_ret <- dist_get(align)
    return(dist_ret)
  } else {
    dist_ret <- dist_ret0
    update_rows <- mcmapply(align_update, min_id = split_test_reduced_min_id, name = split_test_reduced_names, MoreArgs = list(dist_mat = dist_ret), SIMPLIFY = FALSE, mc.cores = cores)
    dist_list <- c(list(dist_ret), update_rows)
    dist_concat <- Reduce(align_bind, dist_list)
    return(dist_concat)
  }
}

umap_process <- function(covid_dist, meta_df) {
  
  acc_names = rownames(covid_dist)
  covid_dist <- dist(covid_dist)
  set.seed(2020)
  covid_umap <- uwot::umap(covid_dist, init = "spectral", metric = "cosine", n_neighbors = 50, min_dist = 0.001, spread = 40, local_connectivity = 10, n_threads = cores*2)
  covid_umap_df <- as.data.frame(covid_umap)
  umap_df_final <- data.frame("Accession" = acc_names, "UMAP_1" = covid_umap_df[,1], "UMAP_2" = covid_umap_df[,2])
  umap_df_final <- merge(umap_df_final, meta_df)
  colnames(umap_df_final) <- c("Accession", "UMAP_1", "UMAP_2", "Region", "Country", "Date")
  return(umap_df_final)
  
}

umap_process_heur <- function(align, fasta, new_dist, new_meta, old_umap) {
  
  umap_0 <- old_umap
  umap_0$Accession <- as.character(umap_0$Accession)
  new_length <- length(readDNAStringSet(fasta))
  acc_names <- new_meta$Accession
  new_acc_names <- acc_names[(length(acc_names) - new_length + 1):length(acc_names)]
  meta_adds <- new_meta[(length(acc_names) - new_length + 1):length(acc_names),]
  new_length_num <- as.list(seq(new_length))
  align_new <- align[(length(align) - new_length + 1):length(align)]
  align_minus <- align[1:(length(align) - new_length)]
  align_new_list <- as.list(align_new)
  
  split_test <- mclapply(align_new, function(x) mclapply(align_minus, function(y) split_dist(x, y), mc.cores = cores), mc.cores = cores)
  split_test_reduced <- mclapply(split_test, function(x) Reduce(cbind, x), mc.cores = cores)
  split_test_reduced_mins <- mclapply(split_test_reduced, min, mc.cores = cores)
  split_test_reduced_min_id <- mclapply(split_test_reduced, which.min, mc.cores = cores)
  split_test_reduced_names <- as.list(names(split_test_reduced_mins))
  
  if (max(unlist(split_test_reduced_mins)) > 1e-4) {
    umap_ret <- umap_process(new_dist, meta_data)
    return(umap_ret)
  } else {
    umap_ret <- umap_0
    umap_update <- function(min_id, number) {
      umap_min_row <- umap_ret[min_id[1],]
      umap_min_row <- umap_min_row[,1:3]
      umap_min_row$Accession[1] <- new_acc_names[number]
      min_row_meta <- meta_adds[number,]
      umap_min_final <- merge(umap_min_row, min_row_meta)
      colnames(umap_min_final) <- colnames(umap_ret)
      return(umap_min_final)
    }
    umap_updates <- mcmapply(umap_update, min_id = split_test_reduced_min_id, number = new_length_num, mc.cores = cores, SIMPLIFY = FALSE)
    umap_updates <- c(list(umap_ret), umap_updates)
    umap_ret <- Reduce(rbind, umap_updates)
    return(umap_ret)
  }
}
    

# Read args

args <- commandArgs(trailingOnly = TRUE)

fasta_file <- as.character(args[1])
align_file <- as.character(args[2])
dist_file <- as.character(args[3])
meta_file <- as.character(args[4])
umap_file <- as.character(args[5])
cores <- as.numeric(args[6])

# Load files and format 

con_fasta <- readDNAStringSet(fasta_file)

align_fasta <- readDNAStringSet(align_file)

dist_dt <- fread(dist_file, stringsAsFactors = FALSE)
dist_mat <- as.matrix(dist_dt[, -1, with = FALSE])
colnames(dist_mat) <- colnames(dist_dt[, -1])
rownames(dist_mat) <- as.vector(dist_dt[[1]])
rm(dist_dt)
gc()

meta_dt <- fread(meta_file, stringsAsFactors = FALSE)
meta_df <- as.data.frame(meta_dt)
meta_df_sub <- meta_df[,c("gisaid_epi_isl", "region", "country", "date")]
colnames(meta_df_sub) <- c("Accession", "Region", "Country", "Date")

umap_dt <- fread(umap_file, stringsAsFactors = FALSE)
umap_df <- as.data.frame(umap_preloaded)

# Compute updated alignment

align_update <- align_get(con_fasta, align_fasta)

# Compute updated distance matrix

dist_mat_up <- dist_get_heur(align_update, con_fasta, dist_mat)

# Process metadata with SIGNAL output

new_accessions <- rownames(dist_mat_up)[(nrow(meta_df_sub)+1):nrow(dist_mat_up)]
meta_new <- data.frame("Accession" = new_accessions, "Region" = paste("Novel", seq(1, length(new_accessions), 1)), "Geo_Location" = paste("Novel", seq(1, length(new_accessions), 1)), "Datetime" = rep((unclass(Sys.Date()) - unclass(as.Date("2019-12-01", format = "%Y-%m-%d"))), length(new_accessions)))
new_meta <- rbind(meta_df_sub, meta_new)

