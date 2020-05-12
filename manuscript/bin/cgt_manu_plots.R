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

setwd("../..")

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

# UMAP figures

umap_figs <- umap_plotter(pre_umap)

umap_figs[[1]]

ggsave("manuscript/data/figures/CGT_umap_region_all.png", device = "png", width = 10, height = 7)

umap_figs[[2]]

ggsave("manuscript/data/figures/CGT_umap_country_all.png", device = "png", width = 16, height = 7)

umap_figs[[3]]

ggsave("manuscript/data/figures/CGT_umap_collect_all.png", device = "png", width = 10, height = 7)

# MST figures

meta_df <- meta
mst_list <- pre_mst
graph_m1 <- mst_list[[1]]
graph_m2 <- mst_list[[2]]
graph_m3 <- mst_list[[3]]
lay <- mst_list[[4]]

ggnet_1 <- ggnetwork(graph_m1, layout = lay)
p1 <- ggplot(ggnet_1, aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_few() +
  geom_edges(color = "gray") +
  geom_nodes(aes(color = Region), size = 3) +
  scale_color_manual(name = "", values = kev_palette[1:length(unique(meta_df$Region))]) +
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

p1

ggsave("manuscript/data/figures/CGT_mst_region_all.png", device = "png", width = 10, height = 7)

ggnet_2 <- ggnetwork(graph_m2, layout = lay)

p2 <- ggplot(ggnet_2, aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_few() +
  geom_edges(color = "gray") +
  geom_nodes(aes(color = Country), size = 3) +
  scale_color_manual(name = "", values = qual_vector[1:length(unique(meta_df$Geo_Location))]) +
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

p2

ggsave("manuscript/data/figures/CGT_mst_country_all.png", device = "png", width = 16, height = 7)

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

p3

ggsave("manuscript/data/figures/CGT_mst_collect_all.png", device = "png", width = 10, height = 7)
                                                                                  
# SNP figures

snp_figs <- snp_plotter(pre_snp, meta)

snp_figs[1]

ggsave("manuscript/data/figures/CGT_snp_region_all.png", device = "png", width = 10, height = 7)

snp_figs[2]

ggsave("manuscript/data/figures/CGT_snp_country_all.png", device = "png", width = 16, height = 7)

snp_figs[3]

ggsave("manuscript/data/figures/CGT_snp_collect_all.png", device = "png", width = 10, height = 7)




