library(Biostrings)
library(DECIPHER)
library(stringi)
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
qual_vector = unlist(mapply(brewer.pal, qual_palettes$maxcolors, rownames(qual_palettes)))

# Functions

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

umap_process <- function(covid_dist, meta_df) {
  
  acc_names = rownames(covid_dist)
  covid_dist <- dist(covid_dist)
  set.seed(2020)
  covid_umap <- uwot::umap(covid_dist, init = "spectral", metric = "cosine", n_neighbors = 50, min_dist = 0.001, spread = 40, local_connectivity = 10)
  covid_umap_df <- as.data.frame(covid_umap)
  umap_df_final <- data.frame("Accession" = acc_names, "UMAP_1" = covid_umap_df[,1], "UMAP_2" = covid_umap_df[,2])
  umap_df_final <- merge(umap_df_final, meta_df)
  colnames(umap_df_final) <- c("Accession", "UMAP_1", "UMAP_2", "Region", "Country", "Date")
  return(umap_df_final)
  
}

mst_graph <- function(covid_dist, meta_data) {
  
  g <- graph.adjacency(covid_dist, mode = "undirected", weighted = TRUE, diag = FALSE)
  g_mst <- mst(g)
  acc_ordering <- match(meta_data[,1], names(V(g_mst)[[]]))
  meta_ordered <- meta_data[order(acc_ordering), ]
  
  g_mst_1 <- g_mst
  g_mst_2 <- g_mst
  g_mst_3 <- g_mst
  
  meta_colors_1 <- (kev_palette[1:length(unique(meta_ordered[,(ncol(meta_ordered)-2)]))])[factor(meta_ordered[, (ncol(meta_ordered)-2)])]
  meta_colors_2 <- (qual_vector[1:length(unique(meta_ordered[,(ncol(meta_ordered)-1)]))])[factor(meta_ordered[, (ncol(meta_ordered)-1)])]
  color_ramp <- colorRampPalette(c("dodgerblue2", "white", "firebrick1"))(max(meta_ordered[,(ncol(meta_ordered))]))
  meta_colors_3 <- color_ramp[meta_ordered[,(ncol(meta_ordered))]]
  
  V(g_mst_1)$color <- meta_colors_1
  V(g_mst_2)$color <- meta_colors_2
  V(g_mst_3)$color <- meta_colors_3
  
  V(g_mst_1)$Region <- meta_ordered[,(ncol(meta_ordered)-2)]
  V(g_mst_2)$Country <- meta_ordered[,(ncol(meta_ordered)-1)]
  V(g_mst_3)$Date <- meta_ordered[,(ncol(meta_ordered))]
  
  lay <- layout_with_graphopt(g_mst_1, niter = 1000)
  
  g_mst_list <- list(g_mst_1, g_mst_2, g_mst_3, lay)
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
    align_final <- data.frame("Position" = rep(pos_full, length(align_tables[,1])), "Meta" = align_tables[,1], "Allele" = align_tables[,2], "Freq" = align_tables[,3])
    return(align_final)
  }
  
  meta_nums <- list(1, 2, 3)
  all_tables <- lapply(meta_nums, function(x) lapply(pos, function(y) table_get(align_df, y, x)))
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

  plot_list <- list(p1, p2, p3)
  return(plot_list)
  
}

mst_plotter <- function(mst_list, meta_df) {
  
  graph_m1 <- mst_list[[1]]
  graph_m2 <- mst_list[[2]]
  graph_m3 <- mst_list[[3]]
  lay <- mst_list[[4]]
  
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
  
  plot_list <- list(p1, p2, p3)
  return(plot_list)
  
}

snp_plotter <- function(snp_list, meta_df) {
  
  snps_1 <- snp_list[[1]]
  colnames(snps_1) <- c("Position", "Region", "Allele", "Freq")
  
  snps_2 <- snp_list[[2]]
  colnames(snps_2) <- c("Position", "Country", "Allele", "Freq")
  
  snps_3 <- snp_list[[3]]
  colnames(snps_3) <- c("Position", "Date", "Allele", "Freq")
  
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
    theme(strip.text = element_text(size = 14, face = "bold")) 
  
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
    theme(strip.text = element_text(size = 14, face = "bold"))
  
  int_text = paste("Functionality currently not supported")
  
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
    theme(strip.text = element_text(size = 14, face = "bold")) 

  plot_list <- list(p1, p2, p3)
  return(plot_list)
  
}
