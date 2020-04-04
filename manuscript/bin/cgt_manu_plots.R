library(Biostrings)
library(DECIPHER)
library(stringi)
library(stringr)
library(ape)
library(igraph)
library(reshape2)
library(dplyr)
library(uwot)
library(ggplot2)
library(ggthemes)
library(Cairo)

# Source all functions

source("global.R")

# Load palettes

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

# Load aligned GISAID profile and metadata

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

pre_aligned <- loadRData("data/dec_aligned_filtered_2020-03-31.RData")
pre_meta <- loadRData("data/covid_filtered_meta_2020-03-31.RData")
pre_dist <- loadRData("data/dec_fasta_dist_2020-03-31.RData")
meta <- pre_meta[,c("Accession", "Region", "Geo_Location", "Datetime")]
pre_umap <- loadRData("data/umap_preloaded.RData")
pre_mst <- loadRData("data/mst_preloaded.RData")
pre_snp <- loadRData("data/snps_preloaded.RData")

# UMAP figure

ggplot(data = pre_umap, aes(x = UMAP_1, y = UMAP_2)) +
  theme_few() +
  geom_jitter(aes(fill = Region), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
  scale_fill_manual(name = "", values = kev_palette[1:length(unique(pre_umap$Region))]) +
  geom_point(data = pre_umap[grep("Novel", pre_umap$Region), ], pch = 21, fill = NA, size = 4, colour = "firebrick1", stroke = 4) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.y = element_text(size = 16, face = "bold")) +
  theme(axis.title.x = element_text(size = 16, face = "bold")) +
  theme(legend.title = element_text(size = 15, face = "bold")) +
  theme(legend.text = element_text(size = 14)) +
  theme(aspect.ratio = 0.8)

ggsave("manuscript/data/figures/CGT_umap_region_all.png", device = "png", width = 10, height = 7)
  
  # MST figure

mst_plotter <- function(mst_list, meta_df) {
  
  graph_m1 <- mst_list[[1]]
  graph_m2 <- mst_list[[2]]
  graph_m3 <- mst_list[[3]]
  lay <- mst_list[[4]]
  
  plot.igraph(graph_m1, vertex.label = NA, vertex.size = 4, edge.width = 1, layout = lay, edge.color = "gray25")
  legend("topleft", legend = levels(factor(meta_df$Region)), fill = kev_palette[1:length(unique(meta_df$Region))])
  
  p1 <- recordPlot()
  
  plot.new()
  
  plot.igraph(graph_m2, vertex.label = NA, vertex.size = 4, edge.width = 1, layout = lay, edge.color = "gray25")
  legend("topleft", legend = levels(factor(meta_df$Geo_Location)), fill = qual_vector[1:length(unique(meta_df$Geo_Location))])
  
  p2 <- recordPlot()
  
  plot.new()
  
  plot.igraph(graph_m3, vertex.label = NA, vertex.size = 4, edge.width = 1, layout = lay, edge.color = "gray25")
  legend("topleft", legend = c("Early", "Mid", "Late"), pt.bg = c("dodgerblue2", "white", "firebrick1"), pt.cex = 1, cex = 1, text.font = 1, pch = 21)
  
  p3 <- recordPlot()
  
  plot_list <- list(p1, p2, p3)
  return(plot_list)
  
}

mst_figs <- mst_plotter(pre_mst, meta)

pdf("manuscript/data/figures/CGT_mst_region_all.pdf", width = 10, height = 10)

mst_figs[1]

dev.off()

plot.new()

pdf("manuscript/data/figures/CGT_mst_region_all_legend.pdf", width = 10, height = 10)

plot.new()

legend("topleft", legend = levels(factor(meta$Region)), fill = kev_palette[1:length(unique(meta$Region))], pt.cex = 1, cex = 1, text.font = 2)

dev.off()
                                                                                  
# SNP figure

snp_figs <- snp_plotter(pre_snp, meta)

snp_figs[1]

ggsave("manuscript/data/figures/CGT_snp_region_all.png", device = "png", width = 10, height = 7)

