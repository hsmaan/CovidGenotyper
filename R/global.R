library(Biostrings)
library(DECIPHER)
library(stringr)
library(ape)
library(igraph)
library(reshape2)
library(dplyr)
library(uwot)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(plotly)
library(ggnetwork)
library(parallel)

# Set core usage

cores <- 2

# Load color palettes    

kev_palette <- c(
  "dodgerblue2", "#E31A1C",
  "green4",
  "#6A3D9A", 
  "#FF7F00", 
  "black", "gold1",
  "skyblue2", "#FB9A99", 
  "palegreen2",
  "#CAB2D6", 
  "#FDBF6F", 
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

qual_palettes = brewer.pal.info[brewer.pal.info$category == "qual", ]
gr_colors = colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)]
qual_vector = sample(gr_colors, 200)

# Functions

align_get <- function(fasta, align) {
  
  covid_seq <- readDNAStringSet(fasta)
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
  if (length(covid_seq) > 1) {
    covid_seq <- AlignSeqs(covid_seq, iterations = 0, refinements = 0, processors = NULL)
  }
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
  new_length <- length(readDNAStringSet(fasta))
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
  covid_dist <- as.dist(covid_dist)
  set.seed(2020)
  covid_umap <- uwot::umap(covid_dist, init = "spectral", metric = "cosine", n_neighbors = 50, min_dist = 0.001, spread = 40, local_connectivity = 10, n_threads = cores*2)
  covid_umap_df <- as.data.frame(covid_umap)
  umap_df_final <- data.frame("Accession" = acc_names, "UMAP_1" = covid_umap_df[,1], "UMAP_2" = covid_umap_df[,2])
  umap_df_final <- merge(umap_df_final, meta_df)
  colnames(umap_df_final) <- c("Accession", "UMAP_1", "UMAP_2", "Region", "Country", "Date", "Travel", "Lineage")
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
    umap_ret <- umap_process(new_dist, new_meta)
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

mst_graph <- function(covid_dist, meta_data) {
  
  g <- graph.adjacency(covid_dist, mode = "undirected", weighted = TRUE, diag = FALSE)
  g_mst <- mst(g)
  acc_ordering <- match(meta_data[,1], names(V(g_mst)[[]]))
  meta_ordered <- meta_data[order(acc_ordering), ]
  
  g_mst_1 <- g_mst
  g_mst_2 <- g_mst
  g_mst_3 <- g_mst
  g_mst_4 <- g_mst
  g_mst_5 <- g_mst
  
  meta_colors_1 <- (kev_palette[1:length(unique(meta_ordered[,(ncol(meta_ordered)-4)]))])[factor(meta_ordered[, (ncol(meta_ordered)-4)])]
  meta_colors_2 <- (qual_vector[1:length(unique(meta_ordered[,(ncol(meta_ordered)-3)]))])[factor(meta_ordered[, (ncol(meta_ordered)-3)])]
  color_ramp <- colorRampPalette(c("dodgerblue2", "white", "firebrick1"))(max(meta_ordered[,(ncol(meta_ordered)-2)]))
  meta_colors_3 <- color_ramp[meta_ordered[,(ncol(meta_ordered)-2)]]
  meta_colors_4 <- (qual_vector[1:length(unique(meta_ordered[,(ncol(meta_ordered)-1)]))])[factor(meta_ordered[, (ncol(meta_ordered)-1)])]
  meta_colors_5 <- (qual_vector[1:length(unique(meta_ordered[,(ncol(meta_ordered))]))])[factor(meta_ordered[, (ncol(meta_ordered))])]

  
  V(g_mst_1)$color <- meta_colors_1
  V(g_mst_2)$color <- meta_colors_2
  V(g_mst_3)$color <- meta_colors_3
  V(g_mst_4)$color <- meta_colors_4
  V(g_mst_5)$color <- meta_colors_5
  
  V(g_mst_1)$Region <- meta_ordered[,(ncol(meta_ordered)-4)]
  V(g_mst_2)$Country <- meta_ordered[,(ncol(meta_ordered)-3)]
  V(g_mst_3)$Date <- meta_ordered[,(ncol(meta_ordered)-2)]
  V(g_mst_4)$Travel <- meta_ordered[,(ncol(meta_ordered)-1)]
  V(g_mst_5)$Lineage <- meta_ordered[,(ncol(meta_ordered))]
  
  lay <- layout_with_graphopt(g_mst_1, niter = 1000)
  
  g_mst_list <- list(g_mst_1, g_mst_2, g_mst_3, g_mst_4, g_mst_5, lay)
  return(g_mst_list)
  
}

snps_get <- function(alignment, metadata, positions) {
  
  align <- alignment
  pos <- paste(positions$Pos, positions$Gene, positions$Effect, sep = " ")
  acc <- metadata[,1]
  meta_vars <- metadata[,-1]
  align_df <- as.data.frame(as.matrix(align))
  meta_order <- match(rownames(align_df), acc)
  align_df <- align_df[order(meta_order),]
  align_df$m5 <- meta_vars[,5]
  align_df$m4 <- meta_vars[,4]
  align_df$m3 <- meta_vars[,3]
  align_df$m2 <- meta_vars[,2]
  align_df$m1 <- meta_vars[,1] # Ensure ordering, cbind() likely reorders rownames and will reorder metadata
  
  freq_pct <- function(d_col) {
    freq_table <- table(d_col)
    pct_table <- as.table(sapply(freq_table, function(x) x/sum(freq_table)))
    return(pct_table)
  }
  
  table_get <- function(df, position, meta_num) {
    pos_full <- position
    pos_num <- as.numeric(str_split_fixed(position, " ", 3)[,1][1])
    align_pos <- df[,c(pos_num, (ncol(df) - as.numeric(meta_num) + 1))] # Backward ordering
    colnames(align_pos) <- c("position", "meta")
    align_grouped <- group_by(.data = align_pos, meta) 
    align_tables <- as.data.frame(do(.data = align_grouped, data.frame(val = freq_pct(.[,1]))))
    align_final <- data.frame(
      "Position" = rep(pos_full, length(align_tables[,1])), 
      "Meta" = align_tables[,1], 
      "Allele" = align_tables[,2], 
      "Freq" = align_tables[,3]
    )
    return(align_final)
  }
  
  meta_nums <- list(1, 2, 3, 4, 5)
  all_tables <- mclapply(meta_nums, function(x) lapply(pos, function(y) table_get(align_df, y, x)), mc.cores = cores)
  table_concat <- lapply(all_tables, function(x) base::Reduce(rbind, x))
  return(table_concat)
  
}

umap_plotter <- function(umap_df) {
  
  p1 <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    theme_few() +
    geom_jitter(aes(fill = Region), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
    scale_fill_manual(name = "", values = kev_palette[1:length(unique(umap_df$Region))]) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 16, face = "bold")) +
    theme(legend.title = element_text(size = 15, face = "bold")) +
    theme(legend.text = element_text(size = 14)) 
  
  p2 <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    theme_few() +
    geom_jitter(aes(fill = Country), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
    scale_fill_manual(name = "", values = qual_vector[1:length(unique(umap_df$Country))]) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 16, face = "bold")) +
    theme(legend.title = element_text(size = 15, face = "bold")) +
    theme(legend.text = element_text(size = 14)) 
  
  p3 <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    theme_few() +
    geom_jitter(aes(color = Date), size = 3, position = "jitter", alpha = 1) +
    scale_color_gradient(name = "Days from \nfirst case", low = "#b92b27", high = "#1565C0") +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 16, face = "bold")) +
    theme(legend.title = element_text(size = 15, face = "bold")) +
    theme(legend.text = element_text(size = 14)) 
  
  p4 <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    theme_few() +
    geom_jitter(aes(fill = Travel), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
    scale_fill_manual(name = "", values = qual_vector[1:length(unique(umap_df$Travel))]) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 16, face = "bold")) +
    theme(legend.title = element_text(size = 15, face = "bold")) +
    theme(legend.text = element_text(size = 14))

  p5 <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    theme_few() +
    geom_jitter(aes(fill = Lineage), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
    scale_fill_manual(name = "", values = qual_vector[1:length(unique(umap_df$Lineage))]) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 16, face = "bold")) +
    theme(legend.title = element_text(size = 15, face = "bold")) +
    theme(legend.text = element_text(size = 14))

  plot_list <- list(p1, p2, p3, p4, p5)
  return(plot_list)
  
}

mst_plotter <- function(mst_list, meta_df) {
  
  graph_m1 <- mst_list[[1]]
  graph_m2 <- mst_list[[2]]
  graph_m3 <- mst_list[[3]]
  graph_m4 <- mst_list[[4]]
  graph_m5 <- mst_list[[5]]
  lay <- mst_list[[6]]
  
  ggnet_1 <- ggnetwork(graph_m1, layout = lay)
  p1 <- ggplot(ggnet_1, aes(x = x, y = y, xend = xend, yend = yend)) +
    theme_few() +
    geom_edges(color = "gray") +
    geom_nodes(aes(fill = Region), size = 3) +
    scale_fill_manual(name = "", values = kev_palette[1:length(unique(meta_df$Region))]) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.line.x = element_blank()) +
    theme(axis.line.y = element_blank()) +
    theme(legend.title = element_text(size = 15, face = "bold")) +
    theme(legend.text = element_text(size = 14)) 
  
  ggnet_2 <- ggnetwork(graph_m2, layout = lay)
  p2 <- ggplot(ggnet_2, aes(x = x, y = y, xend = xend, yend = yend)) +
    theme_few() +
    geom_edges(color = "gray") +
    geom_nodes(aes(fill = Country), size = 3) +
    scale_fill_manual(name = "", values = qual_vector[1:length(unique(meta_df$Geo_Location))]) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.line.x = element_blank()) +
    theme(axis.line.y = element_blank()) +
    theme(legend.title = element_text(size = 15, face = "bold")) +
    theme(legend.text = element_text(size = 14)) 
  
  ggnet_3 <- ggnetwork(graph_m3, layout = lay)
  p3 <- ggplot(ggnet_3, aes(x = x, y = y, xend = xend, yend = yend)) +
    theme_few() +
    geom_edges(color = "gray") +
    geom_nodes(aes(color = Date), size = 3) +
    scale_color_continuous(name = "Days from \nfirst case", low = "#b92b27", high = "#1565C0") + 
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.line.x = element_blank()) +
    theme(axis.line.y = element_blank()) +
    theme(legend.title = element_text(size = 15, face = "bold")) +
    theme(legend.text = element_text(size = 14)) 
  
  ggnet_4 <- ggnetwork(graph_m4, layout = lay)
  p4 <- ggplot(ggnet_4, aes(x = x, y = y, xend = xend, yend = yend)) +
    theme_few() +
    geom_edges(color = "gray") +
    geom_nodes(aes(fill = Travel), size = 3) +
    scale_fill_manual(name = "", values = qual_vector[1:length(unique(meta_df$Country_Exposure))]) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.line.x = element_blank()) +
    theme(axis.line.y = element_blank()) +
    theme(legend.title = element_text(size = 15, face = "bold")) +
    theme(legend.text = element_text(size = 14))
  
  ggnet_5 <- ggnetwork(graph_m5, layout = lay)
  p5 <- ggplot(ggnet_5, aes(x = x, y = y, xend = xend, yend = yend)) +
    theme_few() +
    geom_edges(color = "gray") +
    geom_nodes(aes(fill = Lineage), size = 3) +
    scale_fill_manual(name = "", values = qual_vector[1:length(unique(meta_df$Pangolin_Lineage))]) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.line.x = element_blank()) +
    theme(axis.line.y = element_blank()) +
    theme(legend.title = element_text(size = 15, face = "bold")) +
    theme(legend.text = element_text(size = 14))
  
  plot_list <- list(p1, p2, p3, p4, p5)
  return(plot_list)
  
}

snp_plotter <- function(snp_list, meta_df) {
  
  snps_1 <- snp_list[[1]]
  colnames(snps_1) <- c("Position", "Region", "Allele", "Freq")
  
  snps_2 <- snp_list[[2]]
  colnames(snps_2) <- c("Position", "Country", "Allele", "Freq")
  
  snps_3 <- snp_list[[3]]
  colnames(snps_3) <- c("Position", "Date", "Allele", "Freq")
  
  snps_4 <- snp_list[[4]]
  colnames(snps_4) <- c("Position", "Travel", "Allele", "Freq")

  snps_5 <- snp_list[[5]]
  colnames(snps_5) <- c("Position", "Lineage", "Allele", "Freq")
  
  p1 <- ggplot(data = snps_1, aes(x = Allele, y = Freq)) +
    theme_few () +
    geom_bar(stat = "identity", position = "dodge2", aes(fill = Region), color = "black") +
    scale_fill_manual(name = "", values = kev_palette[1:length(unique(meta_df$Region))]) +
    facet_wrap(~Position, scales = "free_x") +
    labs(x = "Allele", y = "Frequency") + 
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 16, face = "bold")) +
    theme(legend.title = element_text(size = 15, face = "bold")) +
    theme(legend.text = element_text(size = 14)) +
    theme(strip.text = element_text(size = 12, face = "bold")) 
  
  p2 <- ggplot(data = snps_2, aes(x = Allele, y = Freq)) +
    theme_few () +
    geom_bar(stat = "identity", position = "dodge2", aes(fill = Country), width = 1) +
    scale_fill_manual(name = "", values = qual_vector[1:length(unique(meta_df$Geo_Location))]) +
    facet_wrap(~Position, scales = "free") +
    labs(x = "Allele", y = "Frequency") + 
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 16, face = "bold")) +
    theme(legend.title = element_text(size = 15, face = "bold")) +
    theme(legend.text = element_text(size = 14)) +
    theme(strip.text = element_text(size = 12, face = "bold"))
  
  p3 <- ggplot(data = snps_3, aes(x = Allele, y = Freq)) +
    theme_few () +
    geom_bar(stat = "identity", position = "dodge2", aes(color = Date)) +
    scale_color_continuous(name = "Days from \nfirst case", low = "#b92b27", high = "#1565C0") +
    facet_wrap(~Position, scales = "free_x") +
    labs(x = "Allele", y = "Frequency") + 
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 16, face = "bold")) +
    theme(legend.title = element_text(size = 15, face = "bold")) +
    theme(legend.text = element_text(size = 14)) +
    theme(strip.text = element_text(size = 12, face = "bold")) 
  
  p4 <- ggplot(data = snps_4, aes(x = Allele, y = Freq)) +
    theme_few () +
    geom_bar(stat = "identity", position = "dodge2", aes(fill = Travel), width = 1) +
    scale_fill_manual(name = "", values = qual_vector[1:length(unique(meta_df$Country_Exposure))]) +
    facet_wrap(~Position, scales = "free") +
    labs(x = "Allele", y = "Frequency") + 
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 16, face = "bold")) +
    theme(legend.title = element_text(size = 15, face = "bold")) +
    theme(legend.text = element_text(size = 14)) +
    theme(strip.text = element_text(size = 12, face = "bold"))

  p5 <- ggplot(data = snps_5, aes(x = Allele, y = Freq)) +
    theme_few () +
    geom_bar(stat = "identity", position = "dodge2", aes(fill = Lineage), width = 1) +
    scale_fill_manual(name = "", values = qual_vector[1:length(unique(meta_df$Pangolin_Lineage))]) +
    facet_wrap(~Position, scales = "free") +
    labs(x = "Allele", y = "Frequency") + 
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 16, face = "bold")) +
    theme(legend.title = element_text(size = 15, face = "bold")) +
    theme(legend.text = element_text(size = 14)) +
    theme(strip.text = element_text(size = 12, face = "bold"))


  plot_list <- list(p1, p2, p3, p4, p5)
  return(plot_list)
  
}
