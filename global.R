library(Biostrings)
library(DECIPHER)
library(stringi)
library(stringr)
library(ape)
library(igraph)
library(reshape2)
library(dplyr)
library(uwot)

align_get <- function(fasta, align) {
  
  covid_seq <- readDNAStringSet(fasta)
  covid_align <- align
  covid_seq <- RemoveGaps(covid_seq, removeGaps = "all", processors = NULL)
  fasta_final <- AlignProfiles(covid_align, covid_seq)
  return(fasta_final)
  
}

dist_get <- function(align) {
  
  dec_dist <- dist.dna(as.DNAbin(align), model = "K80", as.matrix = TRUE, pairwise.deletion = FALSE)
  colnames(dec_dist) <- (str_split_fixed(colnames(dec_dist), fixed("."), 2)[,1])
  rownames(dec_dist) <- (str_split_fixed(rownames(dec_dist), fixed("."), 2)[,1])
  return(dec_dist)
  
}

mst_graph <- function(covid_dist, meta_data, vertex_cols, meta_num) {
  
  meta_num <- as.numeric(meta_num)
  g <- graph.adjacency(covid_dist, mode = "undirected", weighted = TRUE, diag = FALSE)
  g_mst <- mst(g)
  acc_ordering <- match(meta_data[,1], names(V(g_mst)[[]]))
  meta_ordered <- meta_data[order(acc_ordering), ]
  meta_colors <- (vertex_cols[1:length(unique(meta_ordered[,(meta_num+1)]))])[factor(meta_ordered[,(meta_num+1)])]
  V(g_mst)$color <- meta_colors
  return(g_mst)
  
}

snps_get <- function(alignment, metadata, meta_num) {
  
  meta_num <- as.numeric(meta_num)
  align <- alignment
  pos <- c(1059,1190, 3037, 17858, 18060, 23403, 25563, 27046)
  acc <- metadata[,1]
  meta_var <- metadata[,(meta_num+1)]
  align_df <- as.data.frame(as.matrix(align))
  meta_order <- match(rownames(align_df), acc)
  align_df <- align_df[order(meta_order),]
  align_df$meta <- meta_var
  
  freq_pct <- function(d_col) {
    freq_table <- table(d_col)
    pct_table <- as.table(sapply(freq_table, function(x) x/sum(freq_table)))
    return(pct_table)
  }
  
  table_get <- function(df, position) {
    align_pos <- df[,c(position, ncol(df))]
    align_grouped <- group_by(.data = align_pos, meta) 
    align_tables <- as.data.frame(do(.data = align_grouped, data.frame(val = freq_pct(.[,1]))))
    align_final <- data.frame("Position" = rep(position, length(align_tables[,1])), "Meta" = align_tables[,1], "Allele" = align_tables[,2], "Freq" = align_tables[,3])
    return(align_final)
  }
  
  all_tables <- lapply(pos, function(x) table_get(align_df, x))
  table_concat <- base::Reduce(rbind, all_tables)
  return(table_concat)
  
}

umap_process <- function(covid_dist) {
  
  covid_dist <- dist(covid_dist)
  set.seed(2020)
  covid_umap <- uwot::umap(covid_dist, init = "spectral", metric = "cosine", n_neighbors = 50, min_dist = 0.001, spread = 30, local_connectivity = 10)
  covid_umap_df <- as.data.frame(covid_umap)
  return(covid_umap_df)
  
}